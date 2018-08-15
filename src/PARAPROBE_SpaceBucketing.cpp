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


#include "PARAPROBE_SpaceBucketing.h"

spacebucket::spacebucket()
{
	mdat = sqb();
	buckets.clear();
}


spacebucket::~spacebucket()
{
	freecontainer();
}


void spacebucket::initcontainer( const aabb3d & container )
{
	//##MK::improve bin definition strategy
	mdat.width = static_cast<apt_real>(SPACEBUCKETING_BINWIDTH);
	apt_real rdim = 0.f;
	size_t idim = 0;

	//identify spatial partitioning
	rdim = ceil(container.xsz / mdat.width);
	idim = static_cast<size_t>(rdim);
	mdat.nx = (idim != 0) ? idim : 1;
	rdim = ceil(container.ysz / mdat.width);
	idim = static_cast<size_t>(rdim);
	mdat.ny = (idim != 0) ? idim : 1;
	rdim = ceil(container.zsz / mdat.width);
	idim = static_cast<size_t>(rdim);
	mdat.nz = (idim != 0) ? idim : 1;

	mdat.nxy = mdat.nx * mdat.ny;
	mdat.nxyz = mdat.nx * mdat.ny * mdat.nz;

	mdat.box = container;

	//fill with empty pointer
	buckets.reserve( mdat.nxyz );
	for( size_t z = 0; z < mdat.nz; ++z ) {	//storing in implicit x,y,z order
		for( size_t y = 0; y < mdat.ny; ++y ) {
			for( size_t x = 0; x < mdat.nx; ++x ) {
				buckets.push_back( NULL );
			}
		}
	}
}


void spacebucket::freecontainer( void )
{
	size_t nb = buckets.size();
	for(size_t b = 0; b < nb; ++b) {
		delete buckets.at(b);
		buckets.at(b) = NULL;
	}

	//##MK::reset also mdat to when attempting to reutilize data structure

	buckets.clear();
}


void spacebucket::add_atom( const pos p )
{
	//is point in box which the buckets discretize?
	if ( p.x >= mdat.box.xmi && p.x <= mdat.box.xmx &&
			p.y >= mdat.box.ymi && p.y <= mdat.box.ymx &&
				p.z >= mdat.box.zmi && p.z <= mdat.box.zmx ) { //most likely case
		//in which bucket?
		size_t xx = floor((p.x - mdat.box.xmi) / mdat.box.xsz * static_cast<apt_xyz>(mdat.nx));
		size_t yy = floor((p.y - mdat.box.ymi) / mdat.box.ysz * static_cast<apt_xyz>(mdat.ny));
		size_t zz = floor((p.z - mdat.box.zmi) / mdat.box.zsz * static_cast<apt_xyz>(mdat.nz));

		size_t thisone = xx + yy*mdat.nx + zz*mdat.nxy;

		if ( buckets.at(thisone) != NULL ) { //most likely case until all buckets have been touched
			buckets.at(thisone)->push_back( p );
		}
		else {
			try {
				buckets.at(thisone) = new vector<pos>;
			}
			catch (bad_alloc &ompcroak) {
				complaining( MASTER, "Allocation error in add atom to spacebucket"); return;
			}
			buckets.at(thisone)->push_back( p );
		}
	}
}


void spacebucket::range_rball_noclear_nosort( const pos p, apt_xyz r,
		vector<pos> & candidates )
{
	//does neither clear prior processing nor sort the candidates output array
	//which buckets intruded by sphere of radius r at origin p ?
	apt_real xsc = static_cast<apt_real>(mdat.nx) / mdat.box.xsz;
	apt_real ysc = static_cast<apt_real>(mdat.ny) / mdat.box.ysz;
	apt_real zsc = static_cast<apt_real>(mdat.nz) / mdat.box.zsz;

	apt_real outofbox;
	outofbox = p.x - r - mdat.box.xmi;
	size_t sxmi = (outofbox >= 0.0) ? static_cast<size_t>(floor(outofbox*xsc)) : 0;
	outofbox =  mdat.box.xmx - (p.x + r);
	size_t sxmx = (outofbox >= 0.0) ? static_cast<size_t>(floor((p.x + r - mdat.box.xmi)*xsc)) : mdat.nx-1; //exclusive access

	outofbox = p.y - r - mdat.box.ymi;
	size_t symi = (outofbox >= 0.0) ? static_cast<size_t>(floor(outofbox*ysc)) : 0;
	outofbox =  mdat.box.ymx - (p.y + r);
	size_t symx = (outofbox >= 0.0) ? static_cast<size_t>(floor((p.y + r - mdat.box.ymi)*ysc)) : mdat.ny-1; //exclusive access

	outofbox = p.z - r - mdat.box.zmi;
	size_t szmi = (outofbox >= 0.0) ? static_cast<size_t>(floor(outofbox*zsc)) : 0;
	outofbox =  mdat.box.zmx - (p.z + r);
	size_t szmx = (outofbox >= 0.0) ? static_cast<size_t>(floor((p.z + r - mdat.box.zmi)*xsc)) : mdat.nz-1; //exclusive access

	//scan only buckets within [simi,simx]
	for(size_t z = szmi; z <= szmx; ++z) {
		for(size_t y = symi; y <= symx; ++y) {
			for(size_t x = sxmi; x <= sxmx; ++x) {
				size_t thisone = x + y*mdat.nx + z*mdat.nxy;
				if ( buckets.at(thisone) != NULL ) {
					for(auto it = buckets.at(thisone)->begin(); it != buckets.at(thisone)->end(); ++it) {
						if ( (SQR(it->x - p.x) + SQR(it->y - p.y) + SQR(it->z - p.y)) <= SQR(r) ) { //slightly more likely case
							candidates.push_back( *it );
						}
					}
				} //done checking all candidates within this bucket
			}
		}
	}
}


erase_log spacebucket::erase_rball( const p3d p, apt_xyz r )
{
	//deletes all points inside spherical volume of radius r about p
	apt_real xsc = static_cast<apt_real>(mdat.nx) / mdat.box.xsz;
	apt_real ysc = static_cast<apt_real>(mdat.ny) / mdat.box.ysz;
	apt_real zsc = static_cast<apt_real>(mdat.nz) / mdat.box.zsz;

	//##MK::center are in the box so even if r = 0.f p.x at most mdat.box.xmx
	size_t sxmi = ( (p.x-r) > mdat.box.xmi ) ? static_cast<size_t>(floor((p.x-r-mdat.box.xmi)*xsc)) : 0;
	size_t sxmx = ((p.x+r) < mdat.box.xmx) ? static_cast<size_t>(floor((p.x+r-mdat.box.xmi)*xsc)) : mdat.nx-1; //exclusive access

	size_t symi = ( (p.y-r) > mdat.box.ymi ) ? static_cast<size_t>(floor((p.y-r-mdat.box.ymi)*ysc)) : 0;
	size_t symx = ((p.y+r) < mdat.box.ymx) ? static_cast<size_t>(floor((p.y+r-mdat.box.ymi)*ysc)) : mdat.ny-1; //exclusive access

	size_t szmi = ( (p.z-r) > mdat.box.zmi ) ? static_cast<size_t>(floor((p.z-r-mdat.box.zmi)*zsc)) : 0;
	size_t szmx = ((p.z+r) < mdat.box.zmx) ? static_cast<size_t>(floor((p.z+r-mdat.box.zmi)*zsc)) : mdat.nz-1; //exclusive access

	erase_log status = erase_log();

//cout << "sxmi/sxmx\t\t" << sxmi << "\t\t" << sxmx << endl;
//cout << "symi/symx\t\t" << symi << "\t\t" << symx << endl;
//cout << "szmi/szmx\t\t" << szmi << "\t\t" << szmx << endl;

	//scan only buckets within [simi,simx]
	vector<pos> keep;
	for(size_t z = szmi; z <= szmx; ++z) {
		for(size_t y = symi; y <= symx; ++y) {
			for(size_t x = sxmi; x <= sxmx; ++x) {
				size_t thisone = x + y*mdat.nx + z*mdat.nxy;

				if ( buckets.at(thisone) != NULL ) {
					keep.clear();
					for(auto it = buckets.at(thisone)->begin(); it != buckets.at(thisone)->end(); ++it) {
						if ( (SQR(it->x - p.x) + SQR(it->y - p.y) + SQR(it->z - p.z)) <= SQR(r) ) {
							status.ncleared++;
							continue;
						}
						else {
							status.nkept++;
							keep.push_back( *it );
						}
					}
					if ( keep.size() > 0 ) {
						buckets.at(thisone)->clear();
						buckets.at(thisone)->assign( keep.begin(), keep.end() );
					}
				} //done checking all candidates within this bucket
			}
		}
	}

	return status;
}


size_t spacebucket::get_memory_consumption()
{
	size_t bytes = sizeof(sqb);
	for( size_t i = 0; i < buckets.size(); ++i ) {
		bytes += (1+3)*8; //pointer itself and internal pointer of vector
		bytes += buckets.at(i)->size() * sizeof(pos);
	}
	return bytes; //static_cast<double>(bytes) / (1024 * 1024); //MB
}
