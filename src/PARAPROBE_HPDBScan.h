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


#ifndef __PARAPROBE_HPDBSCAN_H__
#define __PARAPROBE_HPDBSCAN_H__

#include "PARAPROBE_TAPSimHdl.h"

#ifndef HPDBSCAN_UTIL_H
#define HPDBSCAN_UTIL_H

#include <set>
#include <numeric>
#include <cstdlib>
#include <unistd.h>
#include <unordered_set>

#define THREE_DIMENSIONS	3

typedef p3dm1	HPDBPoint;


/**
 * Atomic Operations
 */
template <typename T>
void atomicMin(T* address, T value)
{
    T previous = __sync_fetch_and_add(address, 0);

    while (previous > value) {
        if  (__sync_bool_compare_and_swap(address, previous, value))
        {
            break;
        }
        else
        {
            previous = __sync_fetch_and_add(address, 0);
        }
    }
}


std::ostream& operator<<(std::ostream& os, const p3dm1& pp3 );


/*
 * Output
 */
template <typename T>
std::ostream& operator<<(std::ostream& os, const std::set<T>& set)
{
    os << "{";
    std::copy(set.begin(), set.end(), std::ostream_iterator<T>(os, ", "));
    os << "}";
    return os;
}

template <typename T>
std::ostream& operator<<(std::ostream& os, const std::vector<T>& vec)
{
    os << "[";
    std::copy(vec.begin(), vec.end(), std::ostream_iterator<T>(os, ", "));
    os << "]";
    return os;
}

template <typename T, typename U>
std::ostream& operator<<(std::ostream& os, const std::pair<T, U>& pair)
{
    os << "[" << pair.first << " : " << pair.second << "]";
    return os;
}

template <typename T, typename U>
std::ostream& operator<<(std::ostream& os, const std::map<T, U>& map)
{
    for (auto& pair : map)
    {
        os << pair << std::endl;
    }
    return os;
}

#endif // UTIL_H

////////////////////////////////////////////////////////////////////////////////////////////////

#ifndef HPDBSCAN_POINTS_H
#define HPDBSCAN_POINTS_H


#define REORDER true

const ssize_t NOT_VISITED = std::numeric_limits<ssize_t>::max();
const ssize_t NOISE       = std::numeric_limits<ssize_t>::max() - 1;


class Pointz
{
    size_t m_size;

    std::vector<size_t>  m_cells;
    std::vector<ssize_t> m_clusters;
    std::vector<HPDBPoint>   m_points;
    std::vector<size_t>  m_reorder;

public:
    Pointz(size_t size);

    /**
     *  Iteration
     */
    inline std::vector<HPDBPoint>::iterator begin()
    {
        return this->m_points.begin();
    }

    inline std::vector<HPDBPoint>::const_iterator begin() const
    {
        return this->m_points.begin();
    }

    inline std::vector<HPDBPoint>::iterator end()
    {
        return this->m_points.end();
    }

    inline std::vector<HPDBPoint>::const_iterator end() const
    {
        return this->m_points.end();
    }

    /**
     *  Access
     */
    inline const size_t cell(size_t index) const
    {
        return this->m_cells[index];
    }

    inline const ssize_t cluster(const size_t index) const
    {
        return std::abs(this->m_clusters[index]);
    }

    inline bool corePoint(const size_t index) const
    {
        return this->m_clusters[index] < 0;
    }

    inline HPDBPoint& operator[](size_t index)
    {
        return this->m_points[index];
    }

    inline const size_t size() const
    {
        return this->m_size;
    }

    /**
     * Modifiers
     */
    inline void cell(size_t index, size_t number)
    {
        this->m_cells[index] = number;
    }

    inline void cluster(const size_t index, ssize_t value, bool core)
    {
        atomicMin(&this->m_clusters[index], core ? -value : value);
    }

    inline void overrideCluster(const size_t index, ssize_t value, bool core)
    {
        this->m_clusters[index] = core ? -value : value;
    }

    /**
     * Operations
     */
    void sortByCell(size_t maxDigits);
    void reorder(size_t maxDigits);
};

/**
 * Output
 */
std::ostream& operator<<(std::ostream& os, const Pointz& points);

#endif // POINTS_H



////////////////////////////////////////////////////////////////////////////////////////////////////////

#ifndef HPDBSCAN_SPACE_H
#define HPDBSCAN_SPACE_H


typedef std::map< size_t, size_t >                    CellCounter;
typedef std::map< size_t, std::pair<size_t, size_t> > CellIndex;

class Space
{
    CellIndex           m_cellIndex;

    std::vector<size_t> m_cells;
    HPDBPoint				m_maximum;
    HPDBPoint				m_minimum;

    Pointz&             m_points;
    size_t              m_total;

    /**
     * Initialization
     */
    void computeCells(float epsilon);
    bool probeDimensions(float epsilon, aabb3d & dims );
    void computeDimensions(float epsilon, aabb3d & dims );


public:
    Space(Pointz& points, float epsilon);

    /**
     * Access
     */
    inline const CellIndex& cellIndex() const
    {
        return this->m_cellIndex;
    }

    inline const std::vector<size_t>& cells() const
    {
        return this->m_cells;
    }

    inline size_t cells(size_t dimension) const
    {
        return this->m_cells[dimension];
    }

    inline const HPDBPoint& max() const
    {
        return this->m_maximum;
    }

    /*inline float max(size_t dimension) const
    {
        return this->m_maximum[dimension];

    }*/

    inline const HPDBPoint& min() const
    {
        return this->m_minimum;
    }

    /*inline float min(size_t dimension) const
    {
        return this->m_minimum[dimension];
    }*/

    inline size_t total() const
    {
        return this->m_total;
    }

    /**
     * Operations
     */
    std::vector<size_t> getNeighbors(const size_t cellId) const;
    size_t regionQuery(const size_t pointIndex, const size_t cell, const std::vector<size_t>& neighborPoints, const float EPS2, std::vector<size_t>& minPointsArea) const;
};

/**
 * Output
 */
std::ostream& operator<<(std::ostream& os, const Space& space);

#endif // SPACE_H



#ifndef HPDBSCAN_RULES_H
#define	HPDBSCAN_RULES_H


class Rules
{
    std::map<ssize_t, ssize_t> m_rules;

public:
    Rules();

    const std::map<ssize_t, ssize_t>::const_iterator begin() const
    {
        return this->m_rules.begin();
    }

    const std::map<ssize_t, ssize_t>::const_iterator end() const
    {
        return this->m_rules.end();
    }

    ssize_t rule(const size_t index) const;
    bool update(const ssize_t first, const ssize_t second);
};

void merge(Rules& omp_out, Rules& omp_in);
#pragma omp declare reduction(merge: Rules: merge(omp_out, omp_in)) initializer(omp_priv(omp_orig))

#endif	// RULES_H



///////////////////////////////////////////////////////////////////////////////
#ifndef HPDBSCAN_HPDBSCAN_H
#define HPDBSCAN_HPDBSCAN_H

struct dbscanres
{
	//##MK::quick and dirty
	size_t nCore;
	size_t nClustered;
	size_t nNoise;
	size_t nClusterFound;

	size_t Nmin;
	apt_real Dmax;

	dbscanres() : nCore(0), nClustered(0), nNoise(0),
			nClusterFound(0), Nmin(0), Dmax(0.f) {}
	dbscanres(const size_t _nCore, const size_t _nClu, const size_t _nNoi,
			const size_t _nTot, const size_t _Nm, const apt_real _Dm) :
		nCore(_nCore), nClustered(_nClu), nNoise(_nNoi),
		nClusterFound(_nTot), Nmin(_Nm), Dmax(_Dm) {}
};

class HPDBSCAN
{
public:
    Pointz m_points;

    /*
     * Internal Operations
     */
    void applyRules(const Rules& rules);
    Rules localDBSCAN(const Space &space, float epsilon, size_t minPoints);
    Pointz readIons(const std::vector<p3dm1> & ions );
    Pointz readFile(const std::string& filename);
    dbscanres summarize( TypeSpecDescrStat const & tskdetails, bool const * interior, sqb & vgrd,
    		vector<p3dm1> & ioncloud, const unsigned int maxtypeid, const unsigned int tid, const unsigned int rid );

//public:
    HPDBSCAN(const std::vector<p3dm1>& ions);
    HPDBSCAN(const std::string& filename);
    HPDBSCAN(float* data, size_t size);
    void scan(float epsilon, size_t minPoints);
    unsigned int flag_as_clustered( const unsigned int typid, const unsigned int maxtypid );
    void writeFile(const std::string& outpuFilePath);
    ~HPDBSCAN();
};

#endif // HPDBSCAN_H


struct HPDBParms
{
    std::string file;
    std::string out;
    size_t      threads;
    size_t      minPoints;
    float       epsilon;

    HPDBParms() :
        file("data.csv"),
        out("out.csv"),
        threads(omp_get_max_threads()),
        minPoints(10), epsilon(0.01) {}
    HPDBParms(const float _eps, const size_t _minp) :
    	file(""), out(""), threads(omp_get_max_threads()),
		minPoints(_minp), epsilon(_eps)  {}
};


struct cl3d
{
	cl3d() : id(0), n(0) {}
	cl3d(const size_t _id, const size_t _n) : id(_id), n(_n) {}

	void set_boundary_precipitate()
	{
		this->n = -1;
	}
	bool has_boundary_contact() const
	{
		if ( this->n != -1 )
			return false;
		else
			return true;
	}
	p3d center_of_gravity() const
	{

		if ( this->members.size() > 0 ) {
			p3d cgrv = p3d();
			size_t n = 0;
			for( auto it = this->members.begin(); it != this->members.end(); ++it, ++n ) {
				cgrv.x += it->x;
				cgrv.y += it->y;
				cgrv.z += it->z;
			}
			cgrv.x /= static_cast<apt_real>(n);
			cgrv.y /= static_cast<apt_real>(n);
			cgrv.z /= static_cast<apt_real>(n);
			return cgrv;
		}
		else
			return p3d();
	}

	vector<p3dm1> members;
	size_t id;
	size_t n;
};


class clusterHdl
{
	//instance organizing thread-local processing of cluster
public:
	clusterHdl();
	~clusterHdl();

	void boundary_contact_analysis( bool const * inside, sqb & vgrd );
	void report_size_distr( const string whichtarget, const string againstwhich,
			const unsigned int tid, const unsigned int runid );
	void report_size_distr_vtk( const string whichtarget, const string againstwhich,
			const unsigned int tid, const unsigned int runid );

	vector<cl3d> precipitates;
};

#endif
