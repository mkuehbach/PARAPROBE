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


#include "PARAPROBE_Datatypes.h"


void aabb3d::scale()
{
	this->xsz = this->xmx - this->xmi;
	this->ysz = this->ymx - this->ymi;
	this->zsz = this->zmx - this->zmi;
}

void aabb3d::blowup( const apt_real f )
{
	this->xmi -= f;
	this->xmx += f;
	this->ymi -= f;
	this->ymx += f;
	this->zmi -= f;
	this->zmx += f;
	this->scale();
}

apt_real aabb3d::diag()
{
	return sqrt(SQR(this->xmx-this->xmi)+SQR(this->ymx-this->ymi)+SQR(this->zmx-this->zmi));
}

ostream& operator<<(ostream& in, aabb3d const & val)
{
	in << val.xmi << ";" << val.xmx << "--" << val.ymi << ";" << val.ymx << "--" << val.zmi << ";" << val.zmx << "----" << val.xsz << ";" << val.ysz << ";" << val.zsz << endl;
	return in;
}

ostream& operator<<(ostream& in, jobreceipt const & val)
{
	in << "JobID/ProcessID/ThreadID/WallClock = " << val.jobid << ";" << val.rank << ";" << val.thread << "\t\t" << val.wallclock << endl;
	return in;
}


ostream& operator<<(ostream& in, p3d const & val)
{
	in << val.x << ";" << val.y << ";" << val.z << endl;
	return in;
}


ostream& operator<<(ostream& in, p3dm1 const & val)
{
	in << val.x << ";" << val.y << ";" << val.z << "--" << val.m << endl;
	return in;
}


inline apt_real v3d::len() const
{
	return (this->SQR_len > EPSILON) ? sqrt(this->SQR_len) : 0.f;
}

ostream& operator<<(ostream& in, v3d const & val)
{
	in << val.u << ";" << val.v << ";" << val.w << "----" << val.SQR_len << endl;
	return in;
}

ostream& operator<<(ostream& in, nbor const & val)
{
	in << "Distance/MarkValue = " << val.d << "\t\t" << val.m << endl;
	return in;
}

ostream& operator<<(ostream& in, tri3d const & val)
{
	in << val.x1 << ";" << val.y1 << ";" << val.z1 << "\n";
	in << val.x2 << ";" << val.y2 << ";" << val.z2 << "\n";
	in << val.x3 << ";" << val.y3 << ";" << val.z3 << endl;
	return in;
}

ostream& operator<<(ostream& in, triref3d const & val)
{
	in << val.v1 << "\t\t" << val.v2 << "\t\t" << val.v3 << endl;
	return in;
}


size_t sqb::where( const p3dm1 p )
{
	//3d implicit x+y*nx+z*nx*ny
	apt_real ix = floor((p.x - this->box.xmi) / (this->box.xmx - this->box.xmi) * static_cast<apt_real>(this->nx));
	apt_real iy = floor((p.y - this->box.ymi) / (this->box.ymx - this->box.ymi) * static_cast<apt_real>(this->ny));
	apt_real iz = floor((p.z - this->box.zmi) / (this->box.zmx - this->box.zmi) * static_cast<apt_real>(this->nz));
	size_t res = static_cast<size_t>(ix) + static_cast<size_t>(iy) * this->nx + static_cast<size_t>(iz) * this->nxy;
	return res;
}


ostream& operator<<(ostream& in, sqb const & val)
{
	in << "BinningNXYZ = " << val.nx << ";" << val.ny << ";" << val.nz << "\t\t" << val.nxy << ";" << val.nxyz << "\t\t" << endl;
	in << "BinningBoundingBox = " << val.box << endl;
	return in;
}


geomodel::geomodel(const apt_real _crb, const apt_real _crt, const apt_real _chcb,
		const apt_real _chct, const apt_real _a, const size_t _N )
{
	this->crB = _crb;
	this->crT = _crt;
	this->chcapB = _chcb;
	this->chcapT = _chct;
	this->a = _a;
	this->N = _N;

	/*
	apt_real c1 = this->crB;
	apt_real c2 = this->crT;
	apt_real c3 = this->chcapB;
	apt_real c4 = this->chcapT;
	apt_real aa = this->a; //in nanometer

	size_t Nuz = ATOMS_PER_UNITCELL_FCC;
	apt_real alpha = static_cast<apt_real>(N) * CUBE(aa);
	alpha *= 6.f / (static_cast<apt_real>(ATOMS_PER_UNITCELL_FCC) * PI);

	//volume factor total volume is top cap volume + conical frustum - bottom cap concavity
	apt_real beta = c4*SQR(c2)+CUBE(c4) + 2.f*(SQR(c1)+c1*c2+SQR(c2)) - (c3*SQR(c1)-CUBE(c3));
	apt_real HH = pow( (alpha / beta), (1.f/3.f) );
	*/

	apt_real c1 = this->crT;
	apt_real c2 = this->crB;
	apt_real c3 = this->chcapT;
	apt_real c4 = this->chcapB;
	apt_real aa = this->a; //in nanometer

	apt_real alpha = (static_cast<apt_real>(N) / static_cast<apt_real>(ATOMS_PER_UNITCELL_FCC)) * CUBE(aa) * 6.f / PI;
	//total volume is top cap volume + conical frustum - bottom cap concavity cap volume
	apt_real beta = c3*(3.f*SQR(c1)+SQR(c3)) + 2.f*(SQR(c2)+c2*c1+SQR(c1)) - c4*(3.f*SQR(c2)+SQR(c4));
	apt_real HH = pow( (alpha / beta), (1.f/3.f) );

	this->H = HH;
	this->rB = _crb * HH;
	this->rT = _crt * HH;
	this->hcapB = _chcb * HH;
	this->hcapT = _chct * HH;

	this->mybox.xmi = -1.f * this->rB;
	this->mybox.xmx = +1.f * this->rB;
	this->mybox.ymi = -1.f * this->rB;
	this->mybox.ymx = +1.f * this->rB;
	this->mybox.zmi = 0.f;
	this->mybox.zmx = this->H + this->hcapT;
	this->mybox.scale();

	cout << "alpha/beta = " << alpha << ";" << beta << endl;
	cout << "H/rB/rT/hcB/hcT = " << this->H << ";" << this->rB << ";" << this->rT << ";" << this->hcapB << ";" << this->hcapT << endl;
}


bool geomodel::is_inside( const p3d p)
{
	//MK::check if position p is inside tip volume
	apt_real rB = this->rB;		//radii of the conical frustum
	apt_real rT = this->rT;
	apt_real hB = this->hcapB;	//height of the bottom and top cap elements
	apt_real hT = this->hcapT;
	apt_real H = this->H;		//height of the conical frustum

	//MK::tip volume is axis-aligned with a right-handed cartesian coordinate system
	//the tip axis is parallel to z, the base of the convex hull about the tip is centered at
	//(0,0,0) such that the tip axis is the coordinate system z axis

	//check first if point is in z range above bottom cap but below top cap, most likely case
	apt_real rZ = rB - p.z/H * (rB-rT);
	apt_real sqrxy = SQR(p.x)+SQR(p.y);
	if ( p.z > hB && p.z < H ) { //frustum pillar part
		if ( sqrxy <= SQR(rZ) )
			return true;
		else
			return false;
	}
	else { //only two cases left: either in lower or upper half
		if ( p.z <= hB ) { //lower half
			if ( sqrxy > SQR(rZ) ) //outside conical frustum cross-section circle
				return false;
			else { //inside circle but might still be outside if inside the concave excluded section of the bottom cap
				if ( hB > EPSILON ) { //r = (a^2+h^2)/2h so catch degenerated case
					apt_real bottomradius = (SQR(rB)+SQR(hB)) / (2.f*hB);
					apt_real hhat = hB - p.z;
					apt_real SQRahat = 2.f*bottomradius*hhat - SQR(hhat);
					if ( sqrxy <= SQRahat )
						return false;
					else
						return true;
				}
				//height of the cap is negligible therefore the bottom of the frustum has negligible concavity in which case we are inside
				return true;
			}
		}
		else { //in upper half, only inside if in spherical cap atop of conical frustum
			if ( hT > EPSILON ) { //like with bottom case
				apt_real topradius = (SQR(rT)+SQR(hT)) / (2.f*hT);
				apt_real hhat = hT - (p.z - H);
				apt_real SQRahat = 2.f*topradius*hhat - SQR(hhat);
				if ( sqrxy > SQRahat )
					return false;
				else
					return true;
			}
			//height of cap is again negligible but this time p.z is above so we are outside
			return false;
		}
	}
}


bool SortSpeciAscComposition( speci & first, speci & second )
{
	return first.c < second.c;
}


solutemodel::solutemodel()
{
	urng.seed( MT19937SEED );
	urng.discard( MT19937WARMUP );

	//##MK::debug building of composition, must be the same as in the range file
	//https://www.aircraftmaterials.com/data/aluminium/1200.html
	//http://arc.nucapt.northwestern.edu/refbase/files/Royset-2005.pdf

	//pure aluminiumdummy isotopes
	speci Al27 = speci( 1.000, 26.8, 1 );
	composition.push_back( Al27 );

	/*speci Al27 = speci( 0.980, 26.8, 1 );
	speci Fe56 = speci( 0.010, 55.93, 2 );
	speci Si28 = speci( 0.008, 27.97, 3 );
	speci Sc45 = speci( 0.002, 44.955, 4 );
	composition.push_back( Al27 );
	composition.push_back( Fe56 );
	composition.push_back( Si28 );*/

	//sort and renormalize
	sort( composition.begin(), composition.end(), SortSpeciAscComposition );

	//as cdf is a copy we can modify now the concentration values to be a cdf
	apt_real renorm = 0.f;
	for(size_t i = 0; i < composition.size(); i++)
		renorm += composition.at(i).c;

	apt_real sum = 0.f;
	for(size_t i = 0; i < composition.size(); i++) {
		sum += composition.at(i).c;
		composition.at(i).c = sum / renorm;
//cout << "i/c " << i << "\t\t" << cdf.at(i).c << "\t\t" << cdf.at(i).typid << endl;
	}
}

solutemodel::~solutemodel() {}


unsigned int solutemodel::get_random_speci_typid()
{
	//##MK::error checks

	//acceptance rejection sampling approach with composition CDF to pick a random ion's mass2charge
	//key idea sort compositions in ascending order based on composition compute CDF value

	//pick uniform random number on [0,1)
	uniform_real_distribution<apt_real> distribution(0.f,1.f);
	apt_real cdfval = distribution(urng);

	//use in the mentality of rejection sampling to get corresponding speci index
	for(auto it = composition.begin(); it != composition.end(); ++it) {
		if ( cdfval > it->c ) { //larger than bin end?
			continue;
		}
		return it->typid; //##MK::because for synthetic tips we do so far an in-place ranging;
	}

	//use speci index to get mass to charge
cout << "Unexpected control flow " << endl;
	return UNKNOWNTYPE;
}


apt_real solutemodel::get_random_speci_mq()
{
	//##MK::error checks

	//acceptance rejection sampling approach with composition CDF to pick a random ion's mass2charge
	//key idea sort compositions in ascending order based on composition compute CDF value

	//pick uniform random number on [0,1)
	uniform_real_distribution<apt_real> distribution(0.f,1.f);
	apt_real cdfval = distribution(urng);

	//use in the mentality of rejection sampling to get corresponding speci index
	for(auto it = composition.begin(); it != composition.end(); ++it) {
		if ( cdfval > it->c ) { //larger than bin end?
			continue;
		}
		return it->m2q;
	}

	//use speci index to get mass to charge
cout << "Unexpected control flow " << endl;
	return 0.f;
}




unitcellaggr::unitcellaggr()
{
	this->a = 0.f;

	this->umin = 0;
	this->umax = 0;
	this->vmin = 0;
	this->vmax = 0;
	this->wmin = 0;
	this->wmax = 0;

	this->a1 = v3d();
	this->a2 = v3d();
	this->a3 = v3d();

	this->base.clear();
}


unitcellaggr::unitcellaggr(const apt_real _a, const aabb3d unitbox, const unsigned int model )
{
	a = _a;

	//initialize cubic base vectors
	a1 = v3d( _a*1.f, 0.f, 0.f);
	a2 = v3d( 0.f, _a*1.f, 0.f);
	a3 = v3d( 0.f, 0.f, _a*1.f);

	if ( model == AL3SC ) {
		//initialize Al3Sc base atoms
		base.push_back( p3d(0.f, 0.f, 0.f) ); //Sc 8x 1/8 = 1
		base.push_back( p3d(0.5, 0.5, 0.f) ); //Al
		base.push_back( p3d(0.f, 0.5, 0.5) ); //Al
		base.push_back( p3d(0.5, 0.f, 0.5) ); //Al 6x 1/2 = 3 Al --> Al3Sc okay
	}
	else { //FCC
		//initialize fcc base atoms
		base.push_back( p3d(0.f, 0.f, 0.f) ); //Al 8x 1/8 = 1
		base.push_back( p3d(0.5, 0.5, 0.f) ); //Al
		base.push_back( p3d(0.f, 0.5, 0.5) ); //Al
		base.push_back( p3d(0.5, 0.f, 0.5) ); //Al 6x 1/2 = 3 Al --> 4 units per EZ okay
	}

	//unitbox gives min/max dimensions in nanometer that we have to fill construct on positive sectors of \mathcal{R}^3
	umin = static_cast<int>(floor(unitbox.xmi / _a));
	umax = static_cast<int>(ceil(unitbox.xmx / _a));
	vmin = static_cast<int>(floor(unitbox.ymi / _a));
	vmax = static_cast<int>(ceil(unitbox.ymx / _a));
	wmin = static_cast<int>(floor(unitbox.zmi / _a));
	wmax = static_cast<int>(ceil(unitbox.zmx / _a));

	//unitbox is axis-aligned to standard orientation 0.0, 0.0, 0.0 Bunge Euler fcc crystal lattice
}

unitcellaggr::~unitcellaggr(){}


p3d unitcellaggr::get_atom(const size_t b, const int u, const int v, const int w)
{
	apt_real uu = static_cast<apt_real>(u);
	apt_real vv = static_cast<apt_real>(v);
	apt_real ww = static_cast<apt_real>(w);

	//##MK::implicit origin at 0,0,0
	p3d res = p3d(
			(base[b].x + uu)*a1.u + (base[b].y + vv)*a2.u + (base[b].z + ww)*a3.u,
			(base[b].x + uu)*a1.v + (base[b].y + vv)*a2.v + (base[b].z + ww)*a3.v,
			(base[b].x + uu)*a1.w + (base[b].y + vv)*a2.w + (base[b].z + ww)*a3.w  );

	return res;
}


void msphere::get_atoms( vector<pos> & out )
{
	//##MK::unrotated Al3Sc ideal unit cell centered at sphere center
	//https://materials.springer.com/isp/crystallographic/docs/sd_1922024
	speci Al27 = speci( 0.00, 26.8, 1 );
	speci Sc45 = speci( 0.00, 44.955, 4 );
	apt_real Al3ScLatticeConstant = 0.410; //nm

	apt_real OriginShift = 0.5 * Al3ScLatticeConstant; //MK::base unit cell center is defined at sphere center

	p3d c = center;
	apt_real r = radius;

//cout << "center/c.x/c.y/c.z/r\t\t" << c.x << ";" << c.y << ";" << c.z << "\t\t" << r << endl;

	//aabb3d window = aabb3d( c.x-r, c.x+r, c.y-r, c.y+r, c.z-r, c.z+r );
	aabb3d window = aabb3d( -1.f*r, +1.f*r, -1.f*r, +1.f*r, -1.f*r, +1.f*r );

	unitcellaggr al3sclatt = unitcellaggr( Al3ScLatticeConstant, window, AL3SC );

//cout << "latt/a\t\t\t" << al3sclatt.a << endl;
//cout << "latt/cx/umi/umx\t\t" << c.x << "\t\t" << al3sclatt.umin << "\t\t" << al3sclatt.umax << endl;
//cout << "latt/cy/vmi/vmx\t\t" << c.y << "\t\t" << al3sclatt.vmin << "\t\t" << al3sclatt.vmax << endl;
//cout << "latt/cz/wmi/wmx\t\t" << c.z << "\t\t" << al3sclatt.wmin << "\t\t" << al3sclatt.wmax << endl;

	for(size_t b = 0; b < al3sclatt.base.size(); ++b) {
		//unsigned int mark = (b == 0) ? 4 : 1; //##MK::first defined is Scandium rest Al
		apt_real mass2charge = (b == 0) ? Sc45.m2q : Al27.m2q;

		for(int w = al3sclatt.wmin; w <= al3sclatt.wmax; ++w) {
			for(int v = al3sclatt.vmin; v <= al3sclatt.vmax; ++v) {
				for(int u = al3sclatt.umin; u <= al3sclatt.umax; ++u) {

					p3d ap = al3sclatt.get_atom(b, u, v, w); //##MK::applying temporary shift of origin

					if ( (SQR(ap.x-OriginShift)+SQR(ap.y-OriginShift)+SQR(ap.z-OriginShift)) <= SQR(r) ) {
						out.push_back( pos(
								ap.x - OriginShift + c.x,
								ap.y - OriginShift + c.y,
								ap.z - OriginShift + c.z, mass2charge) );
					} //atom is inside defined sphere
				}
			}
		}
	} //next base
}


bool SortRadiiDesc( apt_real & first, apt_real & second )
{
	return first > second;
}

secondphasemodel::secondphasemodel(const geomodel & geom,
			const size_t N, const apt_real rm, const apt_real rvar )
{
	urng.seed( MT19937SEED );
	urng.discard( MT19937WARMUP );

	//##MK::rm and rvar are in nanometer however we require the distribution parameter mu and sigma
	vector<apt_real> radii;
	lognormal_distribution<apt_real> lognrm( rm, rvar ); //##MK::is in nanometer!
	for( size_t i = 0; i < N; ++i) {
		//fixed size radius value
		radii.push_back( rm );
		//distributed radius value
		//radii.push_back( lognrm(urng) );
	}

	//sorting by descending radius improves packing process speed
	sort( radii.begin(), radii.end(), SortRadiiDesc );

//##MK::DEBUG#############for( size_t i = 0; i < radii.size(); ++i)  { cout << radii.at(i) << endl; } cout << endl << endl;

	uniform_real_distribution<apt_real> unifrnd(0.f,1.f);
	size_t itermax = 100000; //at most so many attempts to place a particle
	size_t iter = 0;
	size_t pid = 0;

	while( particles.size() < N && iter < itermax ) { //attempt placement
		//pick a place at random
		p3d rnd = p3d(
				geom.mybox.xmi + (unifrnd(urng) * geom.mybox.xsz),
				geom.mybox.ymi + (unifrnd(urng) * geom.mybox.ysz),
				geom.mybox.zmi + (unifrnd(urng) * geom.mybox.zsz) );

//cout << "particlesSize/N/iter/itermax\t\t" << particles.size() << ";" << N << ";" << iter << ";" << itermax << "\t\t" << rnd.x << ";" << rnd.y << ";" << rnd.z << endl;

		//overlapping with existent spherical cluster ?
		//MK::precipitates can protrude outside mybox
		bool contact = false;
		//MK::test until intrusion or touching particle found if any
		for( auto pt = particles.begin(); pt != particles.end(); ++pt ) {
			//apt_real distnow = sqrt( SQR(rnd.x - pt->center.x) + SQR(rnd.y - pt->center.y) + SQR(rnd.z - pt->center.z) ); //##MK::sqrt can be avoided here work with SQR values instead
			//apt_real distcrt = radii.at(pid) + pt->radius;

			apt_real distnow = SQR(rnd.x - pt->center.x) + SQR(rnd.y - pt->center.y) + SQR(rnd.z - pt->center.z);
			apt_real distcrt = SQR(radii.at(pid) + pt->radius);

//cout << "\t\t\t" << distnow << "\t\t" << distcrt << endl;
			if ( distnow <= distcrt ) { //most likely the more particles are placed the more overlap
				contact = true;
				break; //MK::we can break if we found at least one overlap
			}
		}

		//decision
		if ( contact == true ) { //later more and more likely
			iter++;
		}
		else { //place the particle
			particles.push_back( msphere( rnd, radii.at(pid) ));
//cout << "rnd/pid/iter\t\t" << rnd.x << ";" << rnd.y << ";" << rnd.z << "\t\t" << pid << "\t\t" << iter << endl;
			pid++;
			iter = 0; //reset counter because next particle placed will again be allowed only maximum
		}
	} //fill 3d box geom.mybox with N randomly placed particles of mean size rm and variance rvar
}


void secondphasemodel::reportParticles( const unsigned int simid, const int rank )
{
	//##MK::suboptimal... one file per rank
	string fn = "PARAPROBE.SimID." + to_string(simid) + ".Rank." + to_string(rank) + ".SyntheticTipCluster.csv";

	ofstream csvlog;
	csvlog.open(fn.c_str(), ofstream::out | ofstream::trunc);
	if (csvlog.is_open() == true) {
		//header
		csvlog << "ClusterID;BarycenterX;BarycenterY;BarycenterZ;Radius\n";
		csvlog<< ";nm;nm;nm;nm\n";
		csvlog << "ClusterID;BarycenterX;BarycenterY;BarycenterZ;Radius\n";

		//report
		for (size_t cid = 0; cid < particles.size(); ++cid ) {
			csvlog << cid << ";" << particles.at(cid).center.x << ";" << particles.at(cid).center.y << ";" << particles.at(cid).center.z << ";" << particles.at(cid).radius << "\n";
		}

		csvlog.flush();
		csvlog.close();
	}
	else {
		stopping("Unable to write synthetic tip placed particles");
	}
}


void secondphasemodel::reportParticlesVTK( const unsigned int simid, const int rank )
{
	//MK::write VTK file showing barycenter of all synthesized cluster in reconstructed space
	//includes particles outside actual tip
	string vtk_io_fn = "PARAPROBE.SimID." + to_string(simid) + ".Rank." + to_string(rank) + ".SyntheticTipCluster.vtk";

	ofstream vtk;
	vtk.open(vtk_io_fn.c_str(), ofstream::out | ofstream::trunc);
	if (vtk.is_open() == true) {
		//header
		vtk << "# vtk DataFile Version 2.0\n";
		vtk << "PARAPROBE Synthetic tip cluster placed\n";
		vtk << "ASCII\n";
		vtk << "DATASET POLYDATA\n";
		vtk << "POINTS " << particles.size() << " double\n";
		for( auto it = particles.begin(); it != particles.end(); ++it ) {
			vtk << it->center.x << " " << it->center.y << " " << it->center.z << "\n";
		}
		vtk << "\n";
		vtk << "VERTICES " << particles.size() << " " << 2*particles.size() << "\n";
		for( size_t i = 0; i < particles.size(); ++i ) {
			vtk << 1 << " " << i << "\n";
		}
		vtk << "\n";
		vtk << "POINT_DATA " << particles.size() << "\n";
		vtk << "FIELD FieldData 1\n";
		vtk << "Radius 1 " << particles.size() << " float\n";
		for( auto it = particles.begin(); it != particles.end(); ++it ) {
			vtk << it->radius << "\n";
		}
		vtk.flush();
		vtk.close();
		cout << "VTK file of synthesized particles positions written to file" << endl;
	}
	else {
		cout << "VTK file of synthesized particles positions was not written" << endl;
	}
}


ostream& operator<<(ostream& in, occupancy const & val)
{
	in << "Result of voxelization is as follows" << "\n";
	in << "TotalBins\t\t" << val.ntotal << "\n";
	in << "VacuumBins\t\t" << val.nvacuum << "\n";
	in << "SurfaceBins\t\t" << val.nsurface << "\n";
	in << "InsideBins\t\t" << val.ninside << "\n";
	in << "SurfaceIons\t\t" << val.nions_surface << "\n";
	in << "InsideIons\t\t" << val.nions_inside << "\n";
	in << "InsideVolume\t\t" << val.volume_inside << " (nm^3)\n";
	return in;
}
