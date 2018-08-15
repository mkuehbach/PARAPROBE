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


#include "PARAPROBE_PeriodicTable.h"


inline apt_real mqrange::width()
{
	return (this->hi - this->lo);
}


inline bool mqrange::inrange(const apt_real val)
{
	if ( val < this->lo ) //too low
		return false;
	if ( val > this->hi ) //too high
		return false;

	return true; //within [lo,hi]
}


unsigned int PeriodicTable::mq2ionname( const apt_real mq )
{
	//MK::maps mass-to-charge to ion type according to known parts of the periodic table
	if ( rangefile_loaded == true ) { //maps a mass to charge ratio to an ionname
		for ( size_t i = 0; i < MQ2IonName.size(); ++i ) {
			if ( MQ2IonName.at(i).size() > 0 ) { //if at all known mass to charge values for this ion speci
				//##MK::earliest possible reject is when completely out of range
				//if ( mq < MQ2IonName.at(i).at(0).lo )
				//	continue;
				//if ( mq > MQ2IonName.at(i).at(MQ2IonName.at(i).size()-1).hi )
				//	continue;

				//no early reject, so there are candidates in mass-to-charge but which
				for ( size_t j = 0; j < MQ2IonName.at(i).size(); ++j ) {
					if (MQ2IonName.at(i).at(j).inrange( mq ) == false ) //try to reject as early as possible...
						continue;
					else
						return i;
				}
			}
		}
	}
	//not returned yet
	//complaining( "PeriodicTable no rangefile was loaded, assigning UNKNOWNTYPEs instead");
	return UNKNOWNTYPE;
}

unsigned int PeriodicTable::ion_name2type( const string name )
{
	map<string, unsigned int>::iterator it = IonTypes.find( name );
	if ( it != IonTypes.end() )
		return it->second;
	else {
		//MK::referring to a completely unknown type
		return UINT32MX;
	}
}

string PeriodicTable::ion_type2name( const unsigned int type )
{
	map<unsigned int, string>::iterator it = IonNames.find( type );
	if ( it != IonNames.end() )
		return it->second;
	else {
		return "COMPLETELYUNKNOWNTYPE";
	}
}

unsigned int PeriodicTable::get_maxtypeid()
{
	return MaximumTypID;
}


bool PeriodicTable::read_rangefile( string ascii_io_fn )
{
	ifstream rrngfile;
	string rrngline;
	istringstream line;
	string datapiece;
	unsigned int typid = UNKNOWNTYPE;
	IonTypes["Unknown"] = UNKNOWNTYPE;
	IonNames[typid] = "Unknown";
	typid++; //MK::typ numeral identification starts at UNKNOWNTYPE+1 !
	vector<mqrange> tmp;
	MQ2IonName.push_back( tmp );
	MQ2IonName.at(UNKNOWNTYPE).push_back( mqrange(0.0,0.001) ); //small dummy range

	//in PARAPROBE_Numerics.h we define UNKNOWNTYPE to be 0 resulting in Fortran indexing of ion types

	rrngfile.open( ascii_io_fn.c_str(), ifstream::in );
	if ( rrngfile.is_open() == true ) {
		//first read keyword [Ions]
		getline( rrngfile, rrngline );
		string is = rrngline.c_str();
		//istringstream line( rrngline );

		//string is = line.str();
		//eliminate potential carriage return
		if (!is.empty() && is[is.size()-1] == '\r')
			is.erase(is.size()-1);
		string should = "[Ions]";

		if ( is.compare(should) == 0 && rrngfile.good() == true ) {
			getline( rrngfile, rrngline );
			istringstream line( rrngline );
			getline( line, datapiece, '=');
			getline( line , datapiece, '=');
			unsigned int nspecies = atoi( datapiece.c_str() );
			for ( unsigned int s = 0; s < nspecies; ++s ) {
				if ( rrngfile.good() == true ) {
						getline( rrngfile, rrngline );
					istringstream line( rrngline );
					getline( line, datapiece, '='); //#MKcheck keyword string prefix Ion
					getline( line, datapiece, '=');

					string physname = string( datapiece.c_str() );
					//eliminate potential carriage return
					if (!physname.empty() && physname[physname.size()-1] == '\r')
						physname.erase(physname.size()-1);

					//avoid multiple reinsertion of ion
					map<string, unsigned int>::iterator it = IonTypes.find( physname );
					if ( it == IonTypes.end() ) { //only if entry does not yet exists work is necessary
						IonTypes[physname] = typid;
						IonNames[typid] = physname;
						vector<mqrange> tmp;
						MQ2IonName.push_back( tmp );
						typid++;
					}
				}
				else {
					stopping( "Rangefile was not completely read in ion definitions");
					rrngfile.close(); rangefile_loaded = false; return false;
				}
			} //read next speci
		} //done reading ions
		else {
			stopping( "Rangefile was not completely read in ion section");
			rrngfile.close(); rangefile_loaded = false; return false;
		}

		if ( rrngfile.good() == true ) {
			getline( rrngfile, rrngline );
			is = rrngline.c_str();
			//eliminate potential carriage return
			if (!is.empty() && is[is.size()-1] == '\r')
				is.erase(is.size()-1);
			should = "[Ranges]";

			if ( is.compare(should) == 0 ) {
				getline( rrngfile, rrngline );
				istringstream line( rrngline );
				getline( line, datapiece, '=');
				getline( line, datapiece, '=');
				unsigned int nranges = atoi( datapiece.c_str() );
				for ( size_t r = 0; r < nranges; ++r ) {
					if ( rrngfile.good() == true ) {
						getline ( rrngfile, rrngline );
						istringstream line( rrngline );
						getline( line, datapiece, '=');

						getline( line, datapiece, ' ');
#ifdef EMPLOY_SINGLEPRECISION
						apt_real low = stof( datapiece.c_str() ); //##MK::check consistence lo < hi etc...
#else
						apt_real low = stod( datapiece.c_str() );
#endif

						getline( line, datapiece, ' ');
#ifdef EMPLOY_SINGLEPRECISION
						apt_real high = stof( datapiece.c_str() );
#else
						apt_real high = stod( datapiece.c_str() );
#endif

						getline( line, datapiece, ' '); //##MK::skip interpreting volume
						getline( line, datapiece, ':');
						string ionkey = datapiece.c_str();
						if (!ionkey.empty() && ionkey[ionkey.size()-1] == '\r') ionkey.erase(ionkey.size()-1); //eliminate carriage return

						map<string, unsigned int>::iterator it = IonTypes.find( ionkey );
						unsigned int thistypid = it->second;
						MQ2IonName.at(thistypid).push_back( mqrange(low, high) );
					}
					else {
						stopping( "Rangefile was not completely in range definitions");
						rrngfile.close(); rangefile_loaded = false; return false;
					}
				} //next range
			}
		}
		else {
			stopping( "Rangefile was not completely read in ranges section");
			rrngfile.close(); rangefile_loaded = false; return false;
		}

		//sort individual keys
		for ( size_t i = 0; i < MQ2IonName.size(); ++i ) {
			sort( MQ2IonName.at(i).begin(), MQ2IonName.at(i).end(), SortMQRangeAscending );
		}

		cout << "RRNGFILE " << ascii_io_fn << " was loaded successfully with a total of now " << typid << " ion types" << endl;
		cout << "Ranging against the following known ion types" << endl;
		for( auto it = IonTypes.begin(); it != IonTypes.end(); it++ ) {
			cout << "\t\t" << it->first << " mapped to TypeID " << it->second << " cross-check name " << IonNames[it->second] << endl;
		}

		//store maximum typid
		MaximumTypID = typid;

		rrngfile.close(); rangefile_loaded = true; return true;
	}
	else {
		complaining( "Unable to load RRNG file " + ascii_io_fn );
		rrngfile.close(); rangefile_loaded = false; return false;
	}
}


TypeSpecDescrStat::TypeSpecDescrStat()
{
	//target = make_pair("UNKNOWNTYPE", UNKNOWNTYPE);
}


TypeSpecDescrStat::~TypeSpecDescrStat()
{
}


/*
bool TypeSpecDescrStat::define_iontask1(const string command, PeriodicTable const & pse )
{
	//is this at all a non-empty string to parse iontype names from?
	if ( command.empty() == true ) {
		return false;
	}

	//is this at all a properly formatted task string?
	//to parse off and is it of the required formatting "RangeType1,RangeType1,RangeType2" e.g. Al,Sc,Zr
	size_t lastcomma = command.rfind(","); //find target type, it is separated by the first comma
	if ( lastcomma == string::npos ) { //if there is not a last one i.e. none then invalid syntax
		return false;
	}
	string lastkw = command.substr(lastcomma+1,command.length()-1);

	istringstream line( command );
	string datapiece;
	getline( line, datapiece, ',');
	auto it = pse.IonTypes.find( datapiece );
	if ( it == pse.IonTypes.end() ) //such name does not exist!
		return false;
	else
		target = make_pair(it->first, it->second);

	//check if frequently desired case all against all
	getline( line, datapiece, ',');
	if ( datapiece.compare("X") == 0 ) { //include all and done
		for (auto jt = pse.IonTypes.begin(); jt != pse.IonTypes.end(); jt++ ) {
			envcandidates.insert( make_pair(jt->first, jt->second) );
		}
		return true;
	}
	else { //okay more complicated case of no all ions to include desired
		auto it = pse.IonTypes.find( datapiece );
		if ( it != pse.IonTypes.end() ) { //include only if it does not yet exist
			auto jt = envcandidates.find( it->first );
			if ( jt == envcandidates.end() ) {
				envcandidates.insert( make_pair(it->first, it->second) );
			} //else do nothing because we have already added the keyword
		} //else cout << "!X case catching attempt to add non-existent type" << endl;
	}

	//process further ion types to include
	while ( 1 )
	{
		getline( line, datapiece, ',' );
		it = pse.IonTypes.find( datapiece );
		if ( it != pse.IonTypes.end() ) { //include only if it does not yet exist
			auto jt = envcandidates.find( it->first );
			if ( jt == envcandidates.end() ) {
				envcandidates.insert( make_pair(it->first, it->second) );
			} //else do nothing because we have already added the keyword
		} //else cout << "while case catching attempt to add non-existent type" << endl;

		if ( datapiece.compare(lastkw) == 0 ) { //terminate if that was the last keyword
			break;
		}
	}
	return true;
}
*/


bool TypeSpecDescrStat::define_iontask2(const string command, PeriodicTable const & pse )
{
	//is this at all a non-empty string to parse iontype names from?
	if ( command.empty() == true ) {
		return false;
	}
	//is there a single comma only in the command string so that the string can specify at all a valid command?
	if ( std::count(command.begin(), command.end(), ',') != 1 ) {
		return false;
	}

	//yes only one comma, so it is a properly formatted task string?
	//CentralType1,CentralTypeN-1-Envtype1,EnvtypeM-1
	//entire task list as follows Al,Al;AlMn,AlMn --->
	//first task Al against Al
	//second task Al or Mn as central ions versus Al or Mn as neighbors

	//split string at the single comma

	istringstream line( command );
	string datapiece;
	getline( line, datapiece, ',');
	//the right part specifies the central ion types
	for (auto jt = pse.IonTypes.begin(); jt != pse.IonTypes.end(); jt++ ) {
		if ( datapiece.find(jt->first) != string::npos )
			trgcandidates.insert( make_pair(jt->first, jt->second) );
	}
	//at least one element as central ion?
	if ( trgcandidates.size() < 1 )
		return false;

	//the left part specifies the neighboring ion types in the environment
	getline( line, datapiece, ',');
	for (auto jt = pse.IonTypes.begin(); jt != pse.IonTypes.end(); jt++ ) {
		if ( datapiece.find(jt->first) != string::npos )
			envcandidates.insert( make_pair(jt->first, jt->second) );
	}
	//at least one element as env ion?
	if ( envcandidates.size() < 1 )
		return false;

	//seems now to be a valid task
	return true;
}




void parse_tasks( const string command, vector<TypeSpecDescrStat> & these, PeriodicTable const & thispse )
{
	//parses a single command string into a vector of individual task descriptives
	if ( command.empty() == true ) { //get out if no command at all
cout << "Parsing no tasks at all" << endl;
		return;
	}

	int numtasks = std::count( command.begin(), command.end(), ';') + 1;
	stringstream parsethis;
	parsethis << command;
	string datapiece;
	for( int tsk = 0; tsk < numtasks; tsk++ ) {
		getline( parsethis, datapiece, ';');

		cout << "Task " << tsk << "____" << datapiece << "____" << endl;

		TypeSpecDescrStat tmp; //speculative accepting of that new task
		these.push_back( tmp );
		TypeSpecDescrStat & thistask = these.back();
		if ( thistask.define_iontask2( datapiece, thispse ) == true ) { //a valid task
			/* ##MK::deprecated
			*cout << endl;
			*cout << "\t\tCentral ion type--->" << thistask.target.first << "/" << thistask.target.second << endl;
			*/
			for ( auto it = thistask.trgcandidates.begin(); it != thistask.trgcandidates.end(); it++) {
				cout << "\t\tCentral ion type--->" << it->first << "/" << it->second << endl;
			}
			for ( auto it = thistask.envcandidates.begin(); it != thistask.envcandidates.end(); it++) {
				cout << "\t\tEnvironment type--->" << it->first << "/" << it->second << endl;
			}
		}
		else { //an invalid task
			cout << "Kicking an invalid task" << endl;
			these.pop_back();
		}
	}
}
