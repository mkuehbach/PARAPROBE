/*
	Copyright Max-Planck-Institut f\"ur Eisenforschung, GmbH, D\"sseldorf
	Data structure, code design, parallel implementation:
	Markus K\"uhbach, 2017-2018

	Third-party contributions:
	Andrew Breen - sequential Matlab code snippets for reconstruction and EPOS
	Markus G\"otz et al. - HPDBScan
	Kartik Kukreja - path compressed union/find
	Lester Hedges - AABBTree

	PARAPROBE --- is an MPI/OpenMP/SIMD-parallelized tool for efficient scalable
	processing of Atom Probe Tomography data targeting back-end processing.
	
	This file is part of PARAPROBE.

	PARAPROBE is free software: you can redistribute it and/or modify
	it under the terms of the GNU General Public License as published by
	the Free Software Foundation, either version 3 of the License, or
	(at your option) any later version.

 	PARAPROBE is distributed in the hope that it will be useful,
	but WITHOUT ANY WARRANTY; without even the implied warranty of
	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
	GNU General Public License for more details.

	You should have received a copy of the GNU General Public License
	along with paraprobe.  If not, see <http://www.gnu.org/licenses/>.
*/

#include "PARAPROBE_TAPSimHdl.h"

ostream& operator<<(ostream& in, tapsim_mvec3d const & val)
{
	in << val.x << ";" << val.y << ";" << val.z;
	return in;
}


ostream& operator<<(ostream& in, tapsim_node const & val)
{
	in << val.x << ";" << val.y << ";" << val.z << "---" << val.id;
	return in;
}


ostream& operator<<(ostream& in, tapsim_phaseVector const & val)
{
	in << val.t << "---" << val.px << ";" << val.py << ";" << val.pz << "---" << val.vx << ";" << val.vy << ";" << val.vz << "---" << val.tetIndex;
	return in;
}


tapsimHdl::tapsimHdl()
{
}


tapsimHdl::~tapsimHdl()
{
    size_t nions = trajectories.size();
    for(size_t ion = 0; ion < nions; ++ion) {
        delete trajectories.at(ion);
        trajectories.at(ion) = NULL;
    }
    trajectories.clear();
}


string tapsimHdl::get_filename( const string prefix, const unsigned int id )
{
	if ( id < 10 )
		return prefix + "." + "0000000" + to_string(id);
	else if ( id < 100 )
		return prefix + "." + "000000" + to_string(id);
	else if ( id < 1000 )
		return prefix + "." + "00000" + to_string(id);
	else if ( id < 10000 )
		return prefix + "." + "0000" + to_string(id);
	else if ( id < 100000 )
		return prefix + "." + "000" + to_string(id);
	else if ( id < 1000000 )
		return prefix + "." + "00" + to_string(id);
	else if ( id < 10000000 )
		return prefix + "." + "0" + to_string(id);
	else if ( id < 100000000 )
		return prefix + "." + to_string(id);
	else
		return "error";
}


bool tapsimHdl::read_binary_log( const string fn )
{
	//actually reading from the file and pumping into trajectories

	ifstream mixedfile;
	string mixedline;
	istringstream line;
	string datapiece;

cout << "Attempting to read " << fn << endl;

	string keyword = "BINARY";
	unsigned long maxcnt = static_cast<unsigned long>(numeric_limits<unsigned int>::max());
	mixedfile.open( fn.c_str(), ifstream::in );
	if ( mixedfile.is_open() == true ) {
		while( mixedfile.good() == true ) {
			//parse ascii lines until coming to keyword BINARY
			getline( mixedfile, mixedline );

			if ( keyword.compare (mixedline.substr(0, 6)) != 0 ) { //most likely no keyword found
				continue;
			}
			else { //keyword found parse expected number of phaseVector struct elements
				istringstream line( mixedline );
				getline( line, datapiece, ' ');
				getline( line, datapiece, ' ');
				unsigned long cnt = stoul(datapiece.c_str());
				if ( cnt >= maxcnt ) { //cast unsafe
					mixedfile.close();
					return false;
				}
				else { //cast safe we know that next parse to read block of structs
					//process into trajectories
					vector<tapsim_phaseVector>* wbuf = NULL;
					try {
						wbuf = new vector<tapsim_phaseVector>;
						wbuf->reserve( static_cast<size_t>(cnt) );
						trajectories.push_back( NULL );
						trajectories.back() = wbuf;
					}
					catch (bad_alloc &tapsimcroak) {
						stopping("Unable to allocate memory for reading in trajectory data");
						mixedfile.close();
						return false;
					}
					unsigned int nj = static_cast<unsigned int>(cnt);

					for(unsigned int j = 0; j < nj; ++j) { //##MK::given file structure horribly slow
						tapsim_phaseVector buf = tapsim_phaseVector();
						mixedfile.read( reinterpret_cast<char*>(&buf), sizeof(tapsim_phaseVector) );
//cout << "\t\t" << buf << endl;
						wbuf->push_back( buf );
					}
cout << "Read " << cnt << " path points" << endl;
				} //done processing path points for ion
			}
		} //keep processing until end of file

		mixedfile.close();
		string mess = "Reading file " + fn + " was successful";
		reporting( mess );
		return true;
	}
	else {
		string mess = "Unable to load file " + fn;
		complaining( mess );
		return false;
	}
}


bool tapsimHdl::read_binary_pathdata( const string fprefix, const unsigned int logID_incr, const unsigned int logID_e )
{
	//##MK::TAPSim v1.0b default trajectory_data filename format is 8 digit leading zeros

	//##MK::first time step is like so trajectory_data.00000001
	string thisfn = get_filename( fprefix, 1 );
	if ( read_binary_log( thisfn ) != true )
		return false;

	//##MK::thereafter like so trajectory_data.00010000
	for ( unsigned int logID = logID_incr; logID <= logID_e; logID += logID_incr ) {
		thisfn = get_filename( fprefix, logID );
		if ( read_binary_log( thisfn ) != true )
			return false;
	}

	return true;
}

bool tapsimHdl::read_node_geometry( const string fn )
{
	//actually reading from the file and pumping into trajectories
	ifstream mixedfile;
	string mixedline;
	istringstream line;
	string datapiece;

cout << "Attempting to read " << fn << endl;
cout << "sizeof(tapsim_node) " << sizeof(tapsim_node) << " Bytes" << endl;

	string keyword = "BINARY";
	unsigned long maxcnt = static_cast<unsigned long>(numeric_limits<unsigned int>::max());
	mixedfile.open( fn.c_str(), ifstream::in );
	if ( mixedfile.is_open() == true ) {
		//reader header first
		getline( mixedfile, mixedline );

		unsigned long cnt = 0;
		if ( keyword.compare (mixedline.substr(0, 6)) != 0 ) { //most likely no keyword found
			mixedfile.close();
			string mess = fn + " has incorrect header formatting or is ASCII";
			complaining(mess);
			return false;
		}
		else { //keyword found parse expected number of tapsim_node struct elements
			istringstream line( mixedline );
			getline( line, datapiece, ' ');
			getline( line, datapiece, ' ');
			cnt = stoul(datapiece.c_str());
			if ( cnt >= maxcnt ) { //cast unsafe
				mixedfile.close();
				string mess = "Cast top unsigned unsafe";
				complaining(mess);
				return false;
			}
			//skip rest of file
			getline( line, datapiece, ' ');
			getline( line, datapiece, ' ');
		}

		nodes.reserve(static_cast<size_t>(cnt));

cout << "Number of nodes " << cnt << endl;

		unsigned int nj = static_cast<unsigned int>(cnt);
		unsigned int j = 0;
		while( mixedfile.good() == true && j < nj) {
			tapsim_mvec3d rbuf1 = tapsim_mvec3d();
			short rbuf2 = numeric_limits<short>::max();
			mixedfile.read( reinterpret_cast<char*>(&rbuf1), sizeof(tapsim_mvec3d) );
			mixedfile.read( reinterpret_cast<char*>(&rbuf2), sizeof(short) );

			if ( rbuf2 != 0 ) {
//cout << "\t\t" << rbuf1 << endl;
cout << "\t\t" << rbuf2 << endl;
			}

			nodes.push_back( tapsim_node( rbuf1.x, rbuf1.y, rbuf1.z, rbuf2 ) );
			j++;
		} //done processing all nodes path points for ion

cout << "Read " << cnt << " nodes" << endl;

		mixedfile.close();
		string mess = "Reading file " + fn + " was successful";
		reporting( mess );
		return true;
	}
	else {
		string mess = "Unable to load file " + fn;
		complaining( mess );
		return false;
	}
}



void tapsimHdl::write_vtk_pathdata( const unsigned int ionID_s, const unsigned int ionID_e )
{
	//sanity check
	if ( ionID_s >= trajectories.size() )
		return;
	if ( ionID_e >= trajectories.size() )
		return;

	for(unsigned int ion = ionID_s; ion <= ionID_e; ++ion) {
		if ( trajectories.at(ion) != NULL ) { //exists?
			string vtk_io_fn = "PARAPROBE.TAPSimPath.IonID." + to_string(ion) + ".vtk";

			ofstream vtk;
			vtk.open( vtk_io_fn.c_str(), ofstream::trunc );
	
			if ( vtk.is_open() == true ) {
				vtk << "# vtk DataFile Version 2.0\n"; //header
				vtk << "PARAPROBE TAPSim Trajectory Ion " << ion << "\n";
				vtk << "ASCII\n";
				vtk << "DATASET POLYDATA\n";

				vector<tapsim_phaseVector>* rbuf = trajectories.at(ion);
				float dist = 0.f;
				size_t nj = 0;
				size_t ni = rbuf->size();
				for (size_t i = 0; i < ni; ++i) {
					if (i > 0) {
						dist += sqrt(	SQR(rbuf->at(i).px - rbuf->at(i-1).px) +
										SQR(rbuf->at(i).py - rbuf->at(i-1).py) +
										SQR(rbuf->at(i).pz - rbuf->at(i-1).pz) );
					}

					if ( dist < 100.0e-9 ) { //all points no further away than 10nm
						nj++;
					}
					else {
						break;
					}
				}
	
				vtk << "POINTS " << nj << " double\n";
				for (size_t i = 0; i < nj; ++i) {
					vtk << rbuf->at(i).px << " " << rbuf->at(i).py << " " << rbuf->at(i).pz << "\n";
				}
				vtk << "\n";
				/*
				vtk << "VERTICES " << nj << " " << 2*nj << "\n";
				for ( size_t i = 0; i < nj; ++i ) {
					vtk << 1 << " " << i << "\n";
				}
				*/
				vtk << "LINES " << 0 << " " << nj << "\n"; //MK::to connect the lines to a polygon set the first to ni-1 and the last to ni-1
				for ( size_t i = 0; i < nj; ++i ) {
					vtk << i << "\n";
				}
				//MK::ranged ion type as field data for coloring in Paraview
				vtk << "POINT_DATA " << nj << "\n";
				vtk << "FIELD FieldData 2\n";
				vtk << "Velocity 3 " << nj << " float\n";
				for(size_t i = 0; i < nj; ++i) {
					vtk << rbuf->at(i).vx << " " << rbuf->at(i).vx << " " << rbuf->at(i).vx << "\n";
				}
				vtk << "\n";
				vtk << "TetIndex 1 " << nj << " int\n";
				for(size_t i = 0; i < nj; ++i) {
					vtk << rbuf->at(i).tetIndex << "\n";
				}


				vtk << endl;
				vtk.flush();
				vtk.close();
			}
			cout << "VTK PathInformation written out for ion " << ion << endl;
		}
	}
}

void tapsimHdl::write_vtk_nodes()
{
	//sanity check
	if ( nodes.size() < 1 )
		return;

	string vtk_io_fn = "PARAPROBE.TAPSimNodeGeometry.vtk";

	ofstream vtk;
	vtk.open( vtk_io_fn.c_str(), ofstream::trunc );

	size_t ni = nodes.size();
	if ( vtk.is_open() == true ) {
		vtk << "# vtk DataFile Version 2.0\n"; //header
		vtk << "PARAPROBE TAPSim Node Geometry " << ni << "\n";
		vtk << "ASCII\n";
		vtk << "DATASET POLYDATA\n";

		vtk << "POINTS " << ni << " double\n";
		for (size_t i = 0; i < ni; ++i) {
			//float xx = (nodes.at(i).x > TAPSIM_EPSILON) ? log10(nodes.at(i).x) : log10(TAPSIM_EPSILON);
			//float yy = (nodes.at(i).y > TAPSIM_EPSILON) ? log10(nodes.at(i).y) : log10(TAPSIM_EPSILON);
			//float zz = (nodes.at(i).z > TAPSIM_EPSILON) ? log10(nodes.at(i).z) : log10(TAPSIM_EPSILON);
			//vtk << xx << " " << yy << " " << zz << "\n";
			vtk << nodes.at(i).x << " " << nodes.at(i).y << " " << nodes.at(i).z << "\n";
		}
		vtk << "\n";
		vtk << "VERTICES " << ni << " " << 2*ni << "\n";
		for ( size_t i = 0; i < ni; ++i ) {
			vtk << 1 << " " << i << "\n";
		}
		//MK::node type as field data
		vtk << "POINT_DATA " << ni << "\n";
		vtk << "FIELD FieldData 1\n";
		vtk << "NodeType 1 " << ni << " int\n";
		for(size_t i = 0; i < ni; ++i) {
			vtk << nodes.at(i).id << "\n";
		}
		vtk << endl;
		vtk.flush();
		vtk.close();

		cout << "VTK NodeGeometry written out" << endl;
	}
}


void tapsimHdl::characterize_pathdensity()
{
}
