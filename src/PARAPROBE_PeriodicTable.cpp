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


PeriodicTable::PeriodicTable()
{
	MaximumTypID = 0;
	rangefile_loaded = false; //flag to prevent multiple reloads
	load_periodictable();
}


PeriodicTable::~PeriodicTable()
{
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


void PeriodicTable::add_iontype_single_or_molecular( vector<string> const & in, mqrange const & ival )
{
	//use PeriodicTable of elements to find whether in contains useful info to build an iontype single or molecular
	//element names are exclusive capitalized case sensitive assumed
	string cand = "";
cout << "adding new iontype" << endl;
	for ( auto it = NuclidTable.begin(); it != NuclidTable.end(); ++it ) { //NuclidTable construction was in order strings in ascending order so
		//does any token of in contain the NuclidName (it->first)
		//format of a valid rrng range file line is Range1=44.8 50.0 Vol:0.01661 Sc:1 Color:33FFFF
		//MK::so do not analyze the 0,1,2 and last token ! otherwise Radon and Cobalt will be added as molecular ions =)
		string nuclidkeyword = it->first;
		bool handle_singlechar_keys = ( nuclidkeyword.size() == 1 ) ? true : false;
		for ( auto tk = in.begin()+3; tk != in.end()-1; ++tk) {
			//special case handling necessary for elements with only one character
			//Fe:1 F:1
			string test = *tk;
			if ( handle_singlechar_keys == true ) { //H, C,
				string keyword_in_rrng = nuclidkeyword + ":";
				if ( test.find(keyword_in_rrng) == string::npos ) {
					continue;
				}
				else {
					cand = cand + nuclidkeyword;
cout << "__" << cand << "__" << endl;
					break; //do not inspect further token to avoid multiple counting of ion
				}
			}
			else {
				if ( test.find( nuclidkeyword ) == string::npos ) { //token does not contain a string like "H" or "Yb"
					continue;
				}
				else {
					cand = cand + it->first;
cout << "__" << cand << "__" << endl;
					break; //do not inspect further token to avoid multiple counting of ion
				}
			}
		}
	}

	//check if an IonType with keyword cand exists already
	if ( IonTypes.find( cand ) == IonTypes.end() ) { //no it doesnt so create
cout << "Single/molecular ion type __" << cand << "__ is new and added" << endl;
		unsigned int typid = static_cast<unsigned int>(IonTypes.size());
		//execute after setting typid
		IonTypes[cand] = typid; //automatic incrementing given that UNKNOWNTYPE was added upon PeriodicTable construction!
		IonNames[typid] = cand;
		vector<mqrange> tmp;
		MQ2IonName.push_back( tmp );
		MQ2IonName.at(typid).push_back( ival );
	}
	else { //yes it does, such add range //##MK::if not overlapping
cout << "Single/molecular ion type __" << cand << "__ exists already just adding range" << endl;
		unsigned int typid = IonTypes[cand];
		MQ2IonName.at(typid).push_back( ival);
	}
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
		if ( is.compare(should) == 0 && rrngfile.good() == true ) { //rangefile as it should
			getline( rrngfile, rrngline );
			istringstream line( rrngline );
			getline( line, datapiece, '=');
			getline( line , datapiece, '=');
			unsigned int nspecies = atoi( datapiece.c_str() ); //instead of reading the species one could just skip this "Ions" section
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
			should = "[Ranges]"; //...and work instead with "Ranges" only and use a periodic table to find combinations of single or molecular ions

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


bool PeriodicTable::read_rangefile2( string ascii_io_fn )
{
	//automatic detection of single and molecular ions from a RRNG file
	//we skip the "Ions" section and read instead immediately the ranges these
	ifstream rrngfile;
	string rrngline;
	istringstream line;
	string datapiece;

	//the zeroth/null element is always the unknown type
	IonTypes["Unknown"] = UNKNOWNTYPE;
	IonNames[UNKNOWNTYPE] = "Unknown";
	//MK::typ numeral identification starts at UNKNOWNTYPE+1 !
	vector<mqrange> tmp;
	MQ2IonName.push_back( tmp );
	MQ2IonName.at(UNKNOWNTYPE).push_back( mqrange(0.0,0.001) ); //small dummy range

	//in PARAPROBE_Numerics.h we define UNKNOWNTYPE to be 0 resulting in Fortran indexing of ion types

	rrngfile.open( ascii_io_fn.c_str(), ifstream::in );
	if ( rrngfile.is_open() == true ) {
		//automatic identification of single and molecular ions given the periodic table of elements
		//so first part of rangefile [Ions] can be skipped lines until eof or keyword "[Ranges]" is found whatever next
		string keyword = "[Ranges]";
		while ( rrngfile.good() == true ) {
			//read a line
			getline( rrngfile, rrngline );
			//does it contain the keyword?
			string is = rrngline.c_str();
			if ( is.find(keyword) == string::npos )
				continue; //may cycle up to eof if file is not a proper rrng
			else
				break; //line contains the keyword so break cycling, rrngline remains current
		}

		while ( rrngfile.good() == true ) { //head controlled, so if previous loop reached eof eofbit is set and will not enter
			//if eofbit not yet set will now filter every line that contains relevant ranging pieces of information
			getline( rrngfile, rrngline );
			//split it into its tokens using boost
			vector<string> tokens;
			boost::split(tokens, rrngline, boost::is_any_of(" \t")); //can handle multiple delimiter different than getline
			if ( 	tokens.size() >= 5 &&
					tokens.at(0).find("Range") != string::npos &&
					tokens.at(2).find("Vol") != string::npos ) { //there is range information available
				//get range m/q interval
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
				mqrange currentrange = mqrange(low, high);
				//now we can use a !case sensitive! searching of element string names in the tokens
				add_iontype_single_or_molecular( tokens, currentrange );
			}
			//whether or not there is range information keep on parsing...
		}

		//sort individual ranges for a given key
		for ( size_t i = 0; i < MQ2IonName.size(); ++i ) {
			sort( MQ2IonName.at(i).begin(), MQ2IonName.at(i).end(), SortMQRangeAscending );
		}

		//store maximum typid
		MaximumTypID = IonTypes.size();
		cout << "The maximum typid is " << MaximumTypID << endl;

		cout << "RRNGFILE " << ascii_io_fn << " was loaded successfully with a total of now " << MaximumTypID << " ion types" << "\n";
		cout << "I know the following elements" << "\n";
		for ( auto it = NuclidTable.begin(); it != NuclidTable.end(); ++it ) {
			cout << " " << it->first;
		}
		cout << endl;
		cout << "I have identified the following single/molecular ions from the rrng file" << "\n";
		for( auto it = IonTypes.begin(); it != IonTypes.end(); it++ ) {
			cout << "\t\t" << it->first << "\t\tmapped to TypeID\t\t" << it->second << "\t\tcross-check name\t\t" << IonNames[it->second] << "\n";
		}
		cout << endl;

		rrngfile.close(); rangefile_loaded = true;
		return true;
		//if no rangeinformation was presented no failure but every ion is considered UNKNOWN debug type
	}
	else {
		complaining( "Unable to load RRNG file " + ascii_io_fn );
		rrngfile.close(); rangefile_loaded = false;
		return false;
	}
}



TypeSpecDescrStat::TypeSpecDescrStat()
{
	//target = make_pair("UNKNOWNTYPE", UNKNOWNTYPE);
}


TypeSpecDescrStat::~TypeSpecDescrStat()
{
}





bool TypeSpecDescrStat::define_iontask3(const string command, PeriodicTable const & pse )
{
	//new type identification scheme:
	//Central1,Central2,...,CentralN-NBor1,NBor2,...,NBorM  i.e.
	//minus sign is central-neighbor block separator
	//comma separates individual elements of centrals or targets

	//reject when command is empty?
	if ( command.empty() == true )
		return false;
	//reject when command has no or more than one minus
	//is there a single comma only in the command string so that the string can specify at all a valid command?
	if ( std::count(command.begin(), command.end(), '-') != 1 )
		return false;

	//split at the central nbor
	istringstream line( command );
	string centers;
	getline( line, centers, '-');
	string nbors;
	getline( line, nbors, '-');

	//how many elements of either centers and nbors?
	size_t ncenters = std::count(centers.begin(), centers.end(), ',') + 1; //##MK +1 because if none there is at least potentially one
	istringstream left( centers );
	string datapiece;
	for( size_t i = 0; i < ncenters; ++i ) {
		getline( left, datapiece, ',');
		for (auto jt = pse.IonTypes.begin(); jt != pse.IonTypes.end(); jt++ ) { //check if keyword is exactly one of the existent
			if ( datapiece.compare(jt->first) != 0 )
				continue;
			else
				trgcandidates.insert( make_pair(jt->first, jt->second) );
		}
	}
	//at least one element as central ion?
	if ( trgcandidates.size() < 1 )
		return false;


	size_t nnbors = std::count(nbors.begin(), nbors.end(), ',') + 1;
	istringstream right( nbors );
	for( size_t i = 0; i < nnbors; ++i ) {
		getline( right, datapiece, ',');
		for (auto jt = pse.IonTypes.begin(); jt != pse.IonTypes.end(); jt++ ) { //check if keyword is exactly one of the existent
			if ( datapiece.compare(jt->first) != 0 )
				continue;
			else
				envcandidates.insert( make_pair(jt->first, jt->second) );
		}
	}
	//at least one element as env ion?
	if ( envcandidates.size() < 1 )
		return false;

	//seems now to be a valid task
	return true;
}


bool TypeSpecDescrStat::define_iontask4(const string command, PeriodicTable const & pse )
{
	//new type identification scheme:
	//Central1,Central2,...,CentralN-NBor1,NBor2,...,NBorM  i.e.
	//minus sign is central-neighbor block separator
	//comma separates individual elements of centrals or targets

	//reformulates userdefined type single/molecular ions into standard internal format
	//for instance consider a rrng file with a Range specification "Al:2 O:1 H:1" a common molecular ion
	//the automatic iontype identification during in read_rangefile will identify from this a molecular ion
	//named AlOoHh, now given that the order in which the elements are found from the range file line
	//is defined but dependent on the order in the for loop over the nuclid table the same molecular ion
	//could be named HhAlOo or all permutations
	//the user should not have to worry about this, ie enter AlOH only now throwing this to the identification
	//algorithm will not find O or H because their keys in the NuclidTable are Hh and Oo respectively
	//this is why we have to internally convert the user input into a consistent representation using the
	//internal format

	//reject when command is empty?
	string mess = "";
	if ( command.empty() == true ) {
		mess = "DefineIontask4::Command __" + command + "__is empty";
		reporting( mess );
		return false;
	}
	//reject when command has no or more than one minus
	//is there a single minus separator only in the command string so that the string can specify at all a valid command?
	if ( std::count(command.begin(), command.end(), '-') != 1 ) {
		mess = "DefineIontask4::Command has not only one minus separator __" + command + "__is empty";
		reporting( mess );
		return false;
	}

	//split at the central nbor
	istringstream line( command );
	string centers;
	getline( line, centers, '-');
	string nbors;
	getline( line, nbors, '-');

	//how many elements of either centers and nbors?
	size_t ncenters = std::count(centers.begin(), centers.end(), ',') + 1; //##MK +1 because if none there is at least potentially one
	istringstream left( centers );
	string datapiece;
	vector<string> temp;
	for( size_t i = 0; i < ncenters; ++i ) {
		getline( left, datapiece, ',');
		cout << "___" << datapiece.c_str() << "___" << datapiece.size() << endl;
		string cand = "";
		temp.clear();
		//parsing of single character element name molecular ions requires additional testing logic
		//take for instance Sc, doing only an accept at first find would find sulphur, some thing for He and hydrogen...
		for ( auto it = pse.NuclidTable.begin(); it != pse.NuclidTable.end(); ++it ) {
			//MK::cycling through the NuclidTable in the same order than upon reading the rangefile
			//such molecular ion iontype name identification is consistent
			//cout << "NuclidTable___" << it->first << "___" << datapiece << "___" << datapiece.c_str() << endl;
			if ( datapiece.find( it->first ) == string::npos ) { //datapiece does not contain a string like "Hh" or "Yb"
				continue;
			}
			else {
				temp.push_back( it->first );
			}
		}
		//now that we know all possible elements of the (molecular)ion we need to build a consistent type string
		//and figure out which they actually are example SSc we would have found S, and Sc as temporaries
		//now test if the given molecular ion name "SSc" is exactly only on of the intermediates
		cout << "Temporary elementname candidates for targets are " << endl;
		for( auto kt = temp.begin(); kt != temp.end(); ++kt )
			cout << " __" << *kt << "__";
		cout << endl;
		bool IsItASingleElementIon = false;
		for( auto kt = temp.begin(); kt != temp.end(); ++kt ) { //take for example "SSc" temp will contain "S" and "Sc"
			if ( datapiece.compare( *kt ) != 0 ) { //"S" != "SSc"
				continue;
			}
			else {
				IsItASingleElementIon = true; //consider as a single ion
				cand = cand + *kt;
				break;
			}
		} //neither "Sc" == "SSc" nor "S" == "SSc" so SSc is a molecular ion
		if ( IsItASingleElementIon == false ) { //consider as a molecular ion
			for( auto kt = temp.begin(); kt != temp.end(); ++kt ) {
				cand = cand + *kt;
			}
		}
		auto jt = pse.IonTypes.find( cand );
		if ( jt != pse.IonTypes.end() ) { //molecular ion does not yet exist
			trgcandidates.insert( make_pair(jt->first, jt->second) );
		}
	}
	//at least one element as central ion?
	if ( trgcandidates.size() < 1 ) {
		mess = "DefineIontask4::trgcandidates.size() < 1";
		reporting( mess );
		return false;
	}


	size_t nnbors = std::count(nbors.begin(), nbors.end(), ',') + 1;
	istringstream right( nbors );
	for( size_t i = 0; i < nnbors; ++i ) {
		getline( right, datapiece, ',');
		string cand = "";
		temp.clear();
		for ( auto it = pse.NuclidTable.begin(); it != pse.NuclidTable.end(); ++it ) {
			//MK::cycling through the NuclidTable in the same order than upon reading the rangefile
			//such molecular ion iontype name identification is consistent
			if ( datapiece.find( it->first ) == string::npos ) { //datapiece does not contain a string like "Hh" or "Yb"
				continue;
			}
			else {
				temp.push_back( it->first );
			}
		}
		cout << "Temporary elementname candidates for neighbors are " << endl;
		for( auto kt = temp.begin(); kt != temp.end(); ++kt )
			cout << " __" << *kt << "__";
		cout << endl;
		//single or molecular ion oracle
		bool IsItASingleElementIon = false;
		for( auto kt = temp.begin(); kt != temp.end(); ++kt ) {
			if ( datapiece.compare( *kt ) != 0 ) {
				continue;
			}
			else {
				IsItASingleElementIon = true;
				cand = cand + *kt;
				break;
			}
		} //neither "Sc" == "SSc" nor "S" == "SSc" so SSc is a molecular ion
		if ( IsItASingleElementIon == false ) { //consider as a molecular ion
			for( auto kt = temp.begin(); kt != temp.end(); ++kt ) {
				cand = cand + *kt;
			}
		}

		auto jt = pse.IonTypes.find( cand );
		if ( jt != pse.IonTypes.end() ) { //molecular ion does not yet exist
			envcandidates.insert( make_pair(jt->first, jt->second) );
		}
	}
	//at least one element as env ion?
	if ( envcandidates.size() < 1 ) {
		mess = "DefineIontask4::envcandidates.size() < 1";
		reporting( mess );
		return false;
	}

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
		//if ( thistask.define_iontask3( datapiece, thispse ) == true ) { //a valid task
		if ( thistask.define_iontask4( datapiece, thispse ) == true ) {
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



void PeriodicTable::load_periodictable()
{
	//MK::so far defines the periodic table of elements only
	//is used for reading rangefiles to allow for automatic detection of single and molecular ion species
	//##MK::can be extended to include nuclides and natural abundances
	nuclid dummy = nuclid();
	//Keyword must be two character string the first capital the second low-case
	//this is why elements with a single character name such as "H" hydrogen require a modified
	//key otherwise a string.find on a an iontype string such as "He" would identify
	//"H" and "He" an form HHe which is incorrect, using Hh prevents this and identifies only He for this example
	//correspondingly molecular ions like AlH become AlHh
	NuclidTable.insert( pair<string,nuclid>("H", dummy) );
	NuclidTable.insert( pair<string,nuclid>("He", dummy) );
	NuclidTable.insert( pair<string,nuclid>("Li", dummy) );
	NuclidTable.insert( pair<string,nuclid>("Be", dummy) );
	NuclidTable.insert( pair<string,nuclid>("B", dummy) );
	NuclidTable.insert( pair<string,nuclid>("C", dummy) );
	NuclidTable.insert( pair<string,nuclid>("N", dummy) );
	NuclidTable.insert( pair<string,nuclid>("O", dummy) );
	NuclidTable.insert( pair<string,nuclid>("F", dummy) );
	NuclidTable.insert( pair<string,nuclid>("Ne", dummy) );
	NuclidTable.insert( pair<string,nuclid>("Na", dummy) );
	NuclidTable.insert( pair<string,nuclid>("Mg", dummy) );
	NuclidTable.insert( pair<string,nuclid>("Al", dummy) );
	NuclidTable.insert( pair<string,nuclid>("Si", dummy) );
	NuclidTable.insert( pair<string,nuclid>("P", dummy) );
	NuclidTable.insert( pair<string,nuclid>("S", dummy) );
	NuclidTable.insert( pair<string,nuclid>("Cl", dummy) );
	NuclidTable.insert( pair<string,nuclid>("Ar", dummy) );
	NuclidTable.insert( pair<string,nuclid>("K", dummy) );
	NuclidTable.insert( pair<string,nuclid>("Ca", dummy) );
	NuclidTable.insert( pair<string,nuclid>("Sc", dummy) );
	NuclidTable.insert( pair<string,nuclid>("Ti", dummy) );
	NuclidTable.insert( pair<string,nuclid>("V", dummy) );
	NuclidTable.insert( pair<string,nuclid>("Cr", dummy) );
	NuclidTable.insert( pair<string,nuclid>("Mn", dummy) );
	NuclidTable.insert( pair<string,nuclid>("Fe", dummy) );
	NuclidTable.insert( pair<string,nuclid>("Co", dummy) );
	NuclidTable.insert( pair<string,nuclid>("Ni", dummy) );
	NuclidTable.insert( pair<string,nuclid>("Cu", dummy) );
	NuclidTable.insert( pair<string,nuclid>("Zn", dummy) );
	NuclidTable.insert( pair<string,nuclid>("Ga", dummy) );
	NuclidTable.insert( pair<string,nuclid>("Ge", dummy) );
	NuclidTable.insert( pair<string,nuclid>("As", dummy) );
	NuclidTable.insert( pair<string,nuclid>("Se", dummy) );
	NuclidTable.insert( pair<string,nuclid>("Br", dummy) );
	NuclidTable.insert( pair<string,nuclid>("Kr", dummy) );
	NuclidTable.insert( pair<string,nuclid>("Rb", dummy) );
	NuclidTable.insert( pair<string,nuclid>("Sr", dummy) );
	NuclidTable.insert( pair<string,nuclid>("Y", dummy) );
	NuclidTable.insert( pair<string,nuclid>("Zr", dummy) );
	NuclidTable.insert( pair<string,nuclid>("Nb", dummy) );
	NuclidTable.insert( pair<string,nuclid>("Mo", dummy) );
	NuclidTable.insert( pair<string,nuclid>("Tc", dummy) );
	NuclidTable.insert( pair<string,nuclid>("Ru", dummy) );
	NuclidTable.insert( pair<string,nuclid>("Rh", dummy) );
	NuclidTable.insert( pair<string,nuclid>("Pd", dummy) );
	NuclidTable.insert( pair<string,nuclid>("Ag", dummy) );
	NuclidTable.insert( pair<string,nuclid>("Cd", dummy) );
	NuclidTable.insert( pair<string,nuclid>("In", dummy) );
	NuclidTable.insert( pair<string,nuclid>("Sn", dummy) );
	NuclidTable.insert( pair<string,nuclid>("Sb", dummy) );
	NuclidTable.insert( pair<string,nuclid>("Te", dummy) );
	NuclidTable.insert( pair<string,nuclid>("I", dummy) );
	NuclidTable.insert( pair<string,nuclid>("Xe", dummy) );
	NuclidTable.insert( pair<string,nuclid>("Cs", dummy) );
	NuclidTable.insert( pair<string,nuclid>("Ba", dummy) );
	NuclidTable.insert( pair<string,nuclid>("La", dummy) );
	NuclidTable.insert( pair<string,nuclid>("Ce", dummy) );
	NuclidTable.insert( pair<string,nuclid>("Pr", dummy) );
	NuclidTable.insert( pair<string,nuclid>("Nd", dummy) );
	NuclidTable.insert( pair<string,nuclid>("Pm", dummy) );
	NuclidTable.insert( pair<string,nuclid>("Sm", dummy) );
	NuclidTable.insert( pair<string,nuclid>("Eu", dummy) );
	NuclidTable.insert( pair<string,nuclid>("Gd", dummy) );
	NuclidTable.insert( pair<string,nuclid>("Tb", dummy) );
	NuclidTable.insert( pair<string,nuclid>("Dy", dummy) );
	NuclidTable.insert( pair<string,nuclid>("Ho", dummy) );
	NuclidTable.insert( pair<string,nuclid>("Er", dummy) );
	NuclidTable.insert( pair<string,nuclid>("Tm", dummy) );
	NuclidTable.insert( pair<string,nuclid>("Yb", dummy) );
	NuclidTable.insert( pair<string,nuclid>("Lu", dummy) );
	NuclidTable.insert( pair<string,nuclid>("Hf", dummy) );
	NuclidTable.insert( pair<string,nuclid>("Ta", dummy) );
	NuclidTable.insert( pair<string,nuclid>("W", dummy) );
	NuclidTable.insert( pair<string,nuclid>("Re", dummy) );
	NuclidTable.insert( pair<string,nuclid>("Os", dummy) );
	NuclidTable.insert( pair<string,nuclid>("Ir", dummy) );
	NuclidTable.insert( pair<string,nuclid>("Pt", dummy) );
	NuclidTable.insert( pair<string,nuclid>("Au", dummy) );
	NuclidTable.insert( pair<string,nuclid>("Hg", dummy) );
	NuclidTable.insert( pair<string,nuclid>("Tl", dummy) );
	NuclidTable.insert( pair<string,nuclid>("Pb", dummy) );
	NuclidTable.insert( pair<string,nuclid>("Bi", dummy) );
	NuclidTable.insert( pair<string,nuclid>("Po", dummy) );
	NuclidTable.insert( pair<string,nuclid>("At", dummy) );
	NuclidTable.insert( pair<string,nuclid>("Rn", dummy) );
	NuclidTable.insert( pair<string,nuclid>("Fr", dummy) );
	NuclidTable.insert( pair<string,nuclid>("Ra", dummy) );
	NuclidTable.insert( pair<string,nuclid>("Ac", dummy) );
	NuclidTable.insert( pair<string,nuclid>("Th", dummy) );
	NuclidTable.insert( pair<string,nuclid>("Pa", dummy) );
	NuclidTable.insert( pair<string,nuclid>("U", dummy) );
	NuclidTable.insert( pair<string,nuclid>("Pu", dummy) );
}

