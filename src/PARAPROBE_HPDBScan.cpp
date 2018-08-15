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

#include "PARAPROBE_HPDBScan.h"


//////////////////////////////////////////////////////////////////////////////////////////


const unsigned short DIGITS   = 10; /* digit count, radix sort buckets */
const unsigned int   POWERS[] = {
    1,
    10,
    100,
    1000,
    10000,
    100000,
    1000000,
    10000000,
    100000000,
    1000000000
};

inline const unsigned int power(size_t exp)
{
    return POWERS[exp];
}


Pointz::Pointz(size_t size) :
    m_cells(std::vector<size_t>(size, 0)),
    m_points(std::vector<HPDBPoint>(size, p3dm1())),
    m_size(size)
{
    this->m_clusters.resize(size);
}


void Pointz::sortByCell(size_t maxDigits)
{
    std::cout << "\tSorting Points... " << std::flush;
    if (maxDigits > DIGITS)
    {
        throw std::invalid_argument("epsilon is too small relative to data space.");
    }

    // reset cluster ids
    std::fill(this->m_clusters.begin(), this->m_clusters.end(), std::numeric_limits<ssize_t>::max());

    std::vector<std::vector<size_t> > buckets(maxDigits, std::vector<size_t>(DIGITS));
    std::vector<size_t> cellBuffer(this->size(), 0);
    std::vector<HPDBPoint> pointsBuffer(this->size(), p3dm1());
    std::vector<size_t> reorderBuffer;
    if (REORDER)
    {
        this->m_reorder = std::vector<size_t>(this->size());
        reorderBuffer = std::vector<size_t>(this->size());
        std::iota (this->m_reorder.begin(), this->m_reorder.end(), 0);
    }

    for (size_t j = 0; j < maxDigits; ++j)
    {
        for(size_t f = 0; f < DIGITS; ++f)
        {
            buckets[j][f] = 0;
        }
    }

    #pragma omp parallel for schedule(static)
    for (size_t i = 0; i < this->size(); ++i)
    {
        for (size_t j = 0; j < maxDigits; ++j)
        {
            const size_t base  = power(j);
            const size_t digit = this->m_cells[i] / base % DIGITS;
            #pragma omp atomic
            ++buckets[j][digit];
        }
    }

    #pragma omp parallel for shared(buckets)
    for (size_t j = 0; j < maxDigits; ++j)
    {
        for (size_t f = 1; f < DIGITS; ++f)
        {
            buckets[j][f] += buckets[j][f-1];
        }
    }

    for (size_t j = 0; j < maxDigits; ++j)
    {
        const size_t base = power(j);
        for (size_t i = this->size() - 1; i < std::numeric_limits<size_t>::max(); --i)
        {
            size_t unit  = this->m_cells[i] / base % DIGITS;
            size_t pos = --buckets[j][unit];

            pointsBuffer[pos] = this->m_points[i];

            cellBuffer[pos] = this->m_cells[i];
            if (REORDER)
            {
                reorderBuffer[pos] = this->m_reorder[i];
            }
        }
        std::copy(pointsBuffer.begin(), pointsBuffer.end(), this->m_points.begin());
        std::copy(cellBuffer.begin(), cellBuffer.end(), this->m_cells.begin());
        if (REORDER)
        {
            std::copy(reorderBuffer.begin(), reorderBuffer.end(), this->m_reorder.begin());
        }
    }

    std::cout << "\t[OK]" << std::endl;
}

void Pointz::reorder(size_t maxDigits)
{
    std::cout << "\tSorting Points Back... " << std::flush;
    std::vector<std::vector<size_t> > buckets(maxDigits, std::vector<size_t>(DIGITS));
    std::vector<HPDBPoint> pointsBuffer(this->size(), p3dm1() );
    std::vector<size_t> cellBuffer(this->size(), 0);
    std::vector<size_t> clusterBuffer(this->size(), 0);
    std::vector<size_t> reorderBuffer(this->size(), 0);

    for (size_t j = 0; j < maxDigits; ++j)
    {
        for(size_t f = 0; f < DIGITS; ++f)
        {
            buckets[j][f] = 0;
        }
    }

    #pragma omp parallel for schedule(static)
    for (size_t i = 0; i < this->size(); ++i)
    {
        for (size_t j = 0; j < maxDigits; ++j)
        {
            const size_t base  = power(j);
            const size_t digit = this->m_reorder[i] / base % DIGITS;
            #pragma omp atomic
            ++buckets[j][digit];
        }
    }

    #pragma omp parallel for shared(buckets)
    for (size_t j = 0; j < maxDigits; ++j)
    {
        for (size_t f = 1; f < DIGITS; ++f)
        {
            buckets[j][f] += buckets[j][f-1];
        }
    }

    for (size_t j = 0; j < maxDigits; ++j)
    {
        const size_t base = power(j);
        for (size_t i = this->size() - 1; i < std::numeric_limits<size_t>::max(); --i)
        {
            size_t unit  = this->m_reorder[i] / base % DIGITS;
            size_t pos = --buckets[j][unit];

            pointsBuffer[pos] = this->m_points[i];

            cellBuffer[pos] = this->m_cells[i];
            clusterBuffer[pos] = this->m_clusters[i];
            reorderBuffer[pos] = this->m_reorder[i];
        }
        std::copy(pointsBuffer.begin(), pointsBuffer.end(), this->m_points.begin());
        std::copy(cellBuffer.begin(), cellBuffer.end(), this->m_cells.begin());
        std::copy(clusterBuffer.begin(), clusterBuffer.end(), this->m_clusters.begin());
        std::copy(reorderBuffer.begin(), reorderBuffer.end(), this->m_reorder.begin());
    }

    std::cout << "\t[OK]" << std::endl;
}

std::ostream& operator<<(std::ostream& os, const Pointz& points)
{
    os << "Points: " << std::endl;
    for (const auto& point : points)
    {
        os << "\t" << point << std::endl;
    }
    os << "Clusters: " << std::endl;
    for (size_t index = 0; index < points.size(); ++index)
    {
        os << "\t" << points.cluster(index) << " " << points.corePoint(index) << std::endl;
    }
    return os;
}




///////////////////////////////////////////////////////////////////////////////////////////////////
void mergeCells(CellCounter& omp_in, CellCounter& omp_out)
{
    for(auto pair: omp_in)
    {
        omp_out[pair.first] += pair.second;
    }
}

void vectorMin(HPDBPoint& omp_in, HPDBPoint& omp_out)
{
    omp_out.x = std::min(omp_out.x, omp_in.x);
    omp_out.y = std::min(omp_out.y, omp_in.y);
    omp_out.z = std::min(omp_out.z, omp_in.z);
}


void vectorMax(HPDBPoint& omp_in, HPDBPoint& omp_out)
{
	omp_out.x = std::max(omp_out.x, omp_in.x);
	omp_out.y = std::max(omp_out.y, omp_in.y);
	omp_out.z = std::max(omp_out.z, omp_in.z);
}


#pragma omp declare reduction(mergeCells: CellCounter: mergeCells(omp_in, omp_out)) initializer(omp_priv(CellCounter()))
#pragma omp declare reduction(vectorMax: HPDBPoint: vectorMax(omp_in, omp_out)) initializer(omp_priv(omp_orig))
#pragma omp declare reduction(vectorMin: HPDBPoint: vectorMin(omp_in, omp_out)) initializer(omp_priv(omp_orig))

Space::Space(Pointz &points, float epsilon) :
    m_cells(std::vector<size_t>(THREE_DIMENSIONS, 0)),
    m_maximum(p3dm1(-F32MX,-F32MX,-F32MX,UNKNOWNTYPE)),
    m_minimum(p3dm1(F32MX,F32MX,F32MX,UNKNOWNTYPE)),
    m_points(points),
    m_total(1)
	//##MK::std::vector<float>(THREE_DIMENSIONS, -std::numeric_limits<float>::max())),
	//##MK::std::vector<float>(THREE_DIMENSIONS,  std::numeric_limits<float>::max())),
{
    std::cout << "\tCalculating Cell Space..." << std::endl;

    //##MK::this is a workaround for datasets so large that prod_i=1,3 dims_i/epsilon_i exceeds 10^DIGITS buckets!
    float use_this_eps = epsilon;
    aabb3d dim = aabb3d( -F32MX, F32MX, -F32MX, F32MX, -F32MX, F32MX );

    if ( this->probeDimensions( epsilon, dim) == false )
    	use_this_eps = 0.4; //nm

    this->computeDimensions(use_this_eps, dim);
    this->computeCells(use_this_eps);

    std::cout << "\tm_total " << this->m_total << endl;
    this->m_points.sortByCell(ceil(log10(this->m_total)));
}


void Space::computeCells(float epsilon)
{
    CellCounter cellCounter;

    std::cout << "\tComputing Cells... " << std::flush;
    #pragma omp parallel for reduction(mergeCells: cellCounter)
    for (size_t i = 0; i < this->m_points.size(); ++i)
    {
        size_t cell    = 0;
        size_t cellAcc = 1;


        float minimum, point;
        size_t dim_index;

        minimum = this->m_minimum.x;
        point = this->m_points[i].x;
        dim_index = (size_t) floor((point - minimum) / epsilon);
        cell    += dim_index * cellAcc;
        cellAcc *= this->m_cells[0];

        minimum = this->m_minimum.y;
		point = this->m_points[i].y;
		dim_index = (size_t) floor((point - minimum) / epsilon);
		cell    += dim_index * cellAcc;
		cellAcc *= this->m_cells[1];

		minimum = this->m_minimum.z;
		point = this->m_points[i].z;
		dim_index = (size_t) floor((point - minimum) / epsilon);
		cell    += dim_index * cellAcc;
		cellAcc *= this->m_cells[2];


        this->m_points.cell(i, cell);
        cellCounter[cell] += 1;
    }

    size_t accumulator      = 0;

    for (auto& cell : cellCounter)
    {
        auto& index  = m_cellIndex[cell.first];
        index.first  = accumulator;
        index.second = cell.second;
        accumulator += cell.second;
    }

    std::cout << "\t[OK]" << std::endl;
}


bool Space::probeDimensions(float epsilon, aabb3d & dims )
{
    aabb3d glim = aabb3d( -F32MX, F32MX, -F32MX, F32MX, -F32MX, F32MX );

	//compute dimensions in parallel using threadlocal aabb3d mylim and fuse in critical
	cout << "\tComputing Dimensions... " << endl;

	#pragma omp parallel shared(glim)
    {
    	aabb3d mylim = aabb3d( -F32MX, F32MX, -F32MX, F32MX, -F32MX, F32MX );
    	#pragma omp for nowait
    	for( size_t iter = 0; iter < this->m_points.size(); ++iter ) {
    		const auto& point = this->m_points[iter];

    		mylim.xmi = std::min( mylim.xmi, point.x );
    		mylim.xmx = std::max( mylim.xmx, point.x );
    		mylim.ymi = std::min( mylim.ymi, point.y );
    		mylim.ymx = std::max( mylim.ymx, point.y );
    		mylim.zmi = std::min( mylim.zmi, point.z );
    		mylim.zmx = std::max( mylim.zmx, point.z );
    	}

		#pragma omp critical
    	{
    		glim.xmi = std::min( glim.xmi, mylim.xmi );
    		glim.xmx = std::max( glim.xmx, mylim.xmx );
    		glim.ymi = std::min( glim.ymi, mylim.ymi );
    		glim.ymx = std::max( glim.ymx, mylim.ymx );
    		glim.zmi = std::min( glim.zmi, mylim.zmi );
    		glim.zmx = std::max( glim.zmx, mylim.zmx );
    	}
    }

    dims = glim;

    //compute total number of cells is below
    size_t mtot = 1;
    size_t cells;
    cells = (size_t) ceil( (glim.xmx - glim.xmi) / epsilon );
    mtot *= cells;
    cells = (size_t) ceil( (glim.ymx - glim.ymi) / epsilon);
    mtot *= cells;
    cells = (size_t) ceil( (glim.zmx - glim.zmi) / epsilon);
    mtot *= cells;

    //probe
    if ( ceil(log10(mtot)) > DIGITS )
    	return false;
    else
    	return true;
}


void Space::computeDimensions(float epsilon, aabb3d & dims )
{
    this->m_minimum = p3dm1( dims.xmi, dims.ymi, dims.zmi, UNKNOWNTYPE );
    this->m_maximum = p3dm1( dims.xmx, dims.ymx, dims.zmx, UNKNOWNTYPE );

    size_t cells;
    cells = (size_t) ceil((this->m_maximum.x - this->m_minimum.x) / epsilon);
    this->m_cells[0] = cells;
    this->m_total   *= cells;
    cells = (size_t) ceil((this->m_maximum.y - this->m_minimum.y) / epsilon);
    this->m_cells[1] = cells;
    this->m_total   *= cells;
    cells = (size_t) ceil((this->m_maximum.z - this->m_minimum.z) / epsilon);
    this->m_cells[2] = cells;
    this->m_total   *= cells;

    cout << "\t\tEpsilon " << epsilon << endl;
    cout << "\t\tNX " << this->m_cells[0] << endl;
    cout << "\t\tNY " << this->m_cells[1] << endl;
    cout << "\t\tNZ " << this->m_cells[2] << "  [OK]" << std::endl;
}


std::vector<size_t> Space::getNeighbors(const size_t cellId) const
{
    const CellIndex& cellIdx = this->m_cellIndex;

    std::vector<size_t> neighborCells;
    neighborCells.reserve(std::pow(3, static_cast<size_t>(THREE_DIMENSIONS) ));
    neighborCells.push_back(cellId);

    size_t lowerSpace     = 1;
    size_t currentSpace   = 1;
    size_t numberOfPoints = cellIdx.find(cellId)->second.second;

    // here be dragons!
    size_t current;
    //X
    currentSpace *= this->m_cells[0];
    for (size_t i = 0, end = neighborCells.size(); i < end; ++i)  {
        current = neighborCells[i];
        // check "left" neighbor - a.k.a the cell in the current dimension that has a lower number
        const long int left = current - lowerSpace;
        auto found    = cellIdx.find(left);
        if (current % currentSpace >= lowerSpace) {
            neighborCells.push_back(left);
            numberOfPoints += found != cellIdx.end() ? found->second.second : 0;
        }
        // check "right" neighbor - a.k.a the cell in the current dimension that has a higher number
        const long int right = current + lowerSpace;
        found = cellIdx.find(right);
        if (current % currentSpace < currentSpace - lowerSpace) {
            neighborCells.push_back(right);
            numberOfPoints += found != cellIdx.end() ? found->second.second : 0;
        }
    }
    lowerSpace = currentSpace;
    //Y
    currentSpace *= this->m_cells[1];
    for (size_t i = 0, end = neighborCells.size(); i < end; ++i)  {
        current = neighborCells[i];
        // check "left" neighbor - a.k.a the cell in the current dimension that has a lower number
        const long int left = current - lowerSpace;
        auto found    = cellIdx.find(left);
        if (current % currentSpace >= lowerSpace) {
            neighborCells.push_back(left);
            numberOfPoints += found != cellIdx.end() ? found->second.second : 0;
        }
        // check "right" neighbor - a.k.a the cell in the current dimension that has a higher number
        const long int right = current + lowerSpace;
        found = cellIdx.find(right);
        if (current % currentSpace < currentSpace - lowerSpace) {
            neighborCells.push_back(right);
            numberOfPoints += found != cellIdx.end() ? found->second.second : 0;
        }
    }
    lowerSpace = currentSpace;
    //Z
    currentSpace *= this->m_cells[2];
    for (size_t i = 0, end = neighborCells.size(); i < end; ++i)  {
        current = neighborCells[i];
        // check "left" neighbor - a.k.a the cell in the current dimension that has a lower number
        const long int left = current - lowerSpace;
        auto found    = cellIdx.find(left);
        if (current % currentSpace >= lowerSpace) {
            neighborCells.push_back(left);
            numberOfPoints += found != cellIdx.end() ? found->second.second : 0;
        }
        // check "right" neighbor - a.k.a the cell in the current dimension that has a higher number
        const long int right = current + lowerSpace;
        found = cellIdx.find(right);
        if (current % currentSpace < currentSpace - lowerSpace) {
            neighborCells.push_back(right);
            numberOfPoints += found != cellIdx.end() ? found->second.second : 0;
        }
    }
    lowerSpace = currentSpace;

    std::vector<size_t> neighborPoints;
    neighborPoints.reserve(numberOfPoints);

    for (size_t neighborCell : neighborCells)
    {
        const auto found = cellIdx.find(neighborCell);
        if (found == cellIdx.end())
        {
            continue;
        }

        const std::pair<size_t, size_t>& locator = found->second;
        neighborPoints.resize(neighborPoints.size() + locator.second);
        auto end = neighborPoints.end();
        std::iota(end - locator.second, end, locator.first);
    }

    return neighborPoints;
}

size_t Space::regionQuery(const size_t pointIndex, const size_t cell, const std::vector<size_t>& neighborPoints, const float EPS2, std::vector<size_t>& minPointsArea) const
{
    const HPDBPoint& point = this->m_points[pointIndex];
    size_t clusterId   = pointIndex + 1; // this MUST be a positive number so that atomicMin will result in correct result with set corePoint bit

    for (size_t neighbor: neighborPoints)
    {
        const HPDBPoint& otherPoint = this->m_points[neighbor];
        float offset            = 0.f;
        offset += std::pow(otherPoint.x - point.x, 2);
        offset += std::pow(otherPoint.y - point.y, 2);
        offset += std::pow(otherPoint.z - point.z, 2);

        if (offset <= EPS2)
        {
            minPointsArea.push_back(neighbor);
            size_t neighborCluster = this->m_points.cluster(neighbor);
            if (neighborCluster != NOT_VISITED && this->m_points.corePoint(neighbor))
            {
                clusterId = std::min(clusterId, neighborCluster);
            }
        }
    }

    return clusterId;
}


std::ostream& operator<<(std::ostream& os, const Space& space)
{
    os << "\tSpace: "    << std::endl;
    os << "\tMin:\t"   << space.min().x << ";" << space.min().y << ";" << space.min().z;
    os << "\tMax:\t"   << space.max().x << ";" << space.max().y << ";" << space.max().z;
    os << "\tCell:\t"  << space.cells();
    os << "\tTotal:\t" << space.total() << std::endl;

    return os;
}


////////////////////////////////////////////////////////////////////////////////////////////
void merge(Rules& omp_out, Rules& omp_in)
{
    for (const auto& rule : omp_in)
    {
        omp_out.update(rule.first, rule.second);
    }
}

Rules::Rules()
{
    this->m_rules[NOISE] = 0L;
}

ssize_t Rules::rule(const size_t index) const
{
    const auto& pair = this->m_rules.find(index);
    if (pair != this->m_rules.end())
    {
        return pair->second;
    }

    return std::numeric_limits<ssize_t>::max();
}

bool Rules::update(const ssize_t first, const ssize_t second)
{
    if (first <= second || first >= NOISE)
    {
        return false;
    }

    const auto& pair = this->m_rules.find(first);
    if (pair != this->m_rules.end())
    {
        if (pair->second > second)
        {
            update(pair->second, second);
            this->m_rules[first]  = second;
        }
        else
        {
            update(second, pair->second );
        }
    }
    else
    {
        this->m_rules[first]  = second;
    }

    return true;
}




////////////////////////////////////////////////////////////////////////////////
void HPDBSCAN::applyRules(const Rules& rules)
{
    #pragma omp parallel for
    for (size_t i = 0; i < this->m_points.size(); ++i)
    {
        const bool core    = this->m_points.corePoint(i);
        ssize_t    cluster = this->m_points.cluster(i);
        ssize_t    found   = rules.rule(cluster);

        while (found < NOISE)
        {
            cluster = found;
            found = rules.rule(found);
        }
        this->m_points.overrideCluster(i, cluster, core);
    }
}


Pointz HPDBSCAN::readIons(const std::vector<p3dm1> & ions )
{
	size_t ni = 0;

	//##MK::preprocessing which ions to account for and how many they are in total comes here

	ni = ions.size(); //##MK::DEBUG take all

	//setup buffer of sufficient initial size
	//##MK::try/catch
	Pointz points(ni);

	//transfer ions
	for(size_t i = 0; i < ni; ++i) {
	   points[i] = ions.at(i);

	   //auto& point = points[i];
	   //ss >> buf;	point.x = atof(buf.c_str()); // get the corordinate of the points
	   //ss >> buf;	point.y = atof(buf.c_str());
	   //ss >> buf;	point.z = atof(buf.c_str());
	   //         	point.m = UNKNOWNTYPE;
	}

	std::cout << "\t" << THREE_DIMENSIONS << "\tDimensions " << std::endl;
	std::cout << "\t" << points.size() << "\tPoints " << std::endl;

	return points;
}


Pointz HPDBSCAN::readFile(const std::string& filename)
{
    std::ifstream file(filename);
    if (!file.is_open())
    {
        std::cout << "Could not open file " << filename << std::endl;
        exit(1);
    }

    std::cout << "Reading " << filename << "... " << std::endl;
    std::string line, buf;
    std::stringstream ss;


    // get the first line
    getline(file, line);
    ss.clear();
    ss << line;

    //##MK::hard-coded dimensions for our specific application in APT
    //size_t dimensions = THREE_DIMENSIONS;

    size_t size = 0;
    // get point count
    file.clear();
    file.seekg (0, std::ios::beg);
    while (!file.eof())
    {
        getline(file, line);
        if(!line.length())
        {
            continue;
        }
        ++size;
    }

    size_t i = 0;
    Pointz points(size);

    // read in points
    file.clear();
    file.seekg (0, std::ios::beg);
    while (!file.eof())
    {
        getline(file, line);
        if (!line.length())
        {
            continue;
        }
        ss.clear();
        ss << line;

        auto& point = points[i];
        ss >> buf;	point.x = atof(buf.c_str()); // get the coordinate of the points
        ss >> buf;	point.y = atof(buf.c_str());
        ss >> buf;	point.z = atof(buf.c_str());
                	point.m = UNKNOWNTYPE;

        ++i;
    }
    file.close();

    std::cout << "\t" << THREE_DIMENSIONS << "\tDimensions " << std::endl;
    std::cout << "\t" << points.size() << "\tPoints " << std::endl;

    return points;
}



void HPDBSCAN::writeFile(const std::string& outputPath)
{
    std::cout << "Writing Output File ...";
    std::ofstream outputFile(outputPath);
    //std::ostream_iterator<float> output(outputFile, " ");

    for (size_t i = 0; i < this->m_points.size(); ++i)
    {
    	const auto& p = this->m_points[i];
        outputFile << p.x << ";" << p.y << ";" << p.z << ";" << p.m << "---->" << this->m_points.cluster(i) << std::endl;

        //std::copy(coordinates.begin(), coordinates.end(), output);
        //outputFile << this->m_points.cluster(i) << std::endl;
    }

    outputFile.close();
    std::cout << "[OK]" << std::endl;
}


Rules HPDBSCAN::localDBSCAN(const Space& space, const float epsilon, const size_t minPoints)
{
    const float EPS2 = std::pow(epsilon, 2);
    const CellIndex& cellIdx = space.cellIndex();

    std::cout << "\tScanning Locally... " << std::endl;
    // unpack keys
    std::vector<size_t> keys;
    keys.reserve(cellIdx.size());
    for (auto iter = cellIdx.begin(); iter != cellIdx.end(); ++iter)
    {
        keys.push_back(iter->first);
    }

    Rules rules;

    // local dbscan
    #pragma omp parallel for schedule(dynamic) reduction(merge: rules)
    for (size_t i = 0; i < keys.size(); ++i)
    {
        const auto& cell = keys[i];
        const auto& pointLocator = cellIdx.find(cell)->second;
        const std::vector<size_t> neighborPoints = space.getNeighbors(cell);

        for (size_t point = pointLocator.first; point < pointLocator.first + pointLocator.second; ++point)
        {
            std::vector<size_t> minPointsArea;
            ssize_t clusterId = space.regionQuery(point, cell, neighborPoints, EPS2, minPointsArea);

            if (minPointsArea.size() >= minPoints)
            {
                this->m_points.cluster(point, clusterId, true);

                for (size_t other : minPointsArea)
                {
                    ssize_t otherClusterId = this->m_points.cluster(other);
                    if (this->m_points.corePoint(other))
                    {
                        auto minmax = std::minmax(otherClusterId, clusterId);
                        rules.update(minmax.second, minmax.first);
                    }
                    this->m_points.cluster(other, clusterId, false);
                }
            }
            else if (this->m_points.cluster(point) == NOT_VISITED)
            {
                this->m_points.cluster(point, NOISE, false);
            }
        }
    }

    return rules;
}

dbscanres HPDBSCAN::summarize( TypeSpecDescrStat const & tskdetails, bool* const interior, sqb const vgrd, vector<p3dm1> & ioncloud,
		const unsigned int maxtypeid, const unsigned int tid, const unsigned int rid )
{
    std::unordered_set<size_t> clusters;
    size_t clusterPoints = 0L;
    size_t noisePoints   = 0L;
    size_t corePoints    = 0L;

    for (size_t i = 0; i < this->m_points.size(); ++i)
    {
        size_t cluster = this->m_points.cluster(i);

        if (cluster == 0) { //noise points are not part of a cluster
            ++noisePoints;
        }
        else { //clusterPoints and corePoints define a cluster
        	 clusters.insert(cluster);
            ++clusterPoints;
        }

        if (this->m_points.corePoint(i)) {
            ++corePoints;
        }
    }
    //cout << "\t" << (clusters.size() - (noisePoints > 0)) << "\tClusters" << endl;
    cout << "\t" << clusters.size() << "\tClusters" << endl;
    cout << "\t" << clusterPoints << "\tCluster Points" << endl;
    cout << "\t" << noisePoints << "\tNoise Points" << endl;
    cout << "\t" << corePoints << "\tCore Points" << endl;

    //characterize cluster
    clusterHdl nanogod;
    nanogod.precipitates.reserve( clusters.size() );

    //mapping of HPDBScan cluster ID (which can be nonconsecutive) to index on nanogod.precipitate array
    map<size_t, size_t> cid2idx;
    size_t newlabel = 0;
    for( auto it = clusters.begin(); it != clusters.end(); ++it ) {
    	cid2idx.insert( pair<size_t,size_t>(*it, newlabel) );
    	nanogod.precipitates.push_back( cl3d( newlabel, 0 ) ); //MK::*it gives val of the unordered_set object
//cout << *it << "\t\t" << nanogod.precipitates.size()-1 << "\t\t" << newlabel << endl;
		newlabel++;
    }

    for (size_t i = 0; i < this->m_points.size(); ++i)
    {
        size_t cluster = this->m_points.cluster(i);
        if (cluster == 0) //noise ion
        	continue;
        else { //add ion to its cluster
        	auto it = cid2idx.find(cluster);
        	if ( it != cid2idx.end() ) {
        		nanogod.precipitates.at( it->second ).n++;
        		nanogod.precipitates.at( it->second ).members.push_back( this->m_points[i] );
        		//these.at( cid2idx.at(cluster) ).members.push_back( m_points.at(i) ); //throws
        	}
        }
    }

cout << "Boundary contact analysis of all cluster..." << endl;
    nanogod.boundary_contact_analysis( interior, vgrd );
cout << "Boundary contact analysis of all cluster done" << endl;

    string whichtarg = ""; //tskdetails.target.first;
    for( auto jt = tskdetails.trgcandidates.begin(); jt != tskdetails.trgcandidates.end(); jt++ ) {
    	whichtarg += jt->first;
    }
   	string whichcand = "";
    for( auto jt = tskdetails.envcandidates.begin(); jt != tskdetails.envcandidates.end(); jt++) {
    	whichcand += jt->first;
    }
    nanogod.report_size_distr( whichtarg, whichcand, tid, rid );
    nanogod.report_size_distr_vtk( whichtarg, whichcand, tid, rid );

    //flag ions as clustered
    unsigned int mxtypid = maxtypeid;
    size_t excluded = 0;
    for(size_t i = 0; i < ioncloud.size(); ++i) {
    	size_t clusterid = this->m_points.cluster(i);
    	if ( clusterid != 0 ) { //non-noise ion
    		ioncloud.at(i).m = flag_as_clustered( ioncloud.at(i).m, mxtypid );
    		++excluded;
    	}
    }

    cout << endl << "Flag2ExcludeNoNoise " << excluded << endl;
    //now all ions of icld which sc found to cluster have flagged ion types,
   	//while all others still have their original type


    //wrap up results
    return dbscanres( corePoints, clusterPoints, noisePoints, static_cast<size_t>(clusters.size()), 0, 0.f ); //##MK::last two values dummies
}


HPDBSCAN::HPDBSCAN(const std::vector<p3dm1>& ioncloud) :
	m_points(readIons(ioncloud)) {}


HPDBSCAN::HPDBSCAN(const std::string &filename) :
    m_points(readFile(filename))
{}


inline ssize_t stoll(const char* optarg)
{
    size_t  length;
    ssize_t parsed;

    try
    {
        parsed = std::stoll(optarg, &length);
    }
    catch (const std::logic_error& e)
    {
        parsed = 0L;
    }

    return length < strlen(optarg) ? 0L : parsed;
}

inline float stof(const char* optarg)
{
    size_t length;
    float  parsed;

    try
    {
        parsed = std::stof(optarg, &length);
    }
    catch (const std::logic_error& e)
    {
        parsed = 0.0f;
    }

    return length < strlen(optarg) ? 0.0f : parsed;
}


void HPDBSCAN::scan(float epsilon, size_t minPoints)
{
    Space space(this->m_points, epsilon);
    Rules rules = this->localDBSCAN(space, epsilon, minPoints);
    this->applyRules(rules);
    if (REORDER)
    {
        this->m_points.reorder(ceil(log10(this->m_points.size())));
    }
	//return (this->summarize());
}

unsigned int HPDBSCAN::flag_as_clustered( const unsigned int typid, const unsigned int maxtypid )
{
	if ( typid <= maxtypid )
		return (typid + ION_IN_CLUSTER);
	else
		return typid;
}

HPDBSCAN::~HPDBSCAN()
{
}



clusterHdl::clusterHdl()
{
}

clusterHdl::~clusterHdl()
{
}


void clusterHdl::boundary_contact_analysis( bool* const inside, sqb const & smartgrd )
{
	//MK::go through all my precipitates and check whether they with respect to the binning grd
	//and the binarized 3d indicator inside protrude outside the inside
	for(auto it = precipitates.begin(); it != precipitates.end(); it++) {
		//it->n = it->members.size();

		//for all points forming the cluster check for protrusion outside the inside
		vector<p3dm1> const & these = it->members;

		/*
		for(size_t mb = 0; mb < these.size(); mb++) {
			p3dm1 pp = these.at(mb);
			size_t cxyz = smartgrd.where( pp );
			if ( inside[cxyz] == true ) {
				continue;
			}
			else {
				it->set_boundary_precipitate();
				break;
			}
		} //test members until all done or one protrude outside
		*/

	} //test all precipitates
}


void clusterHdl::report_size_distr( const string whichtarget, const string againstwhich, const unsigned int tid, const unsigned int runid )
{
	string fn = "PARAPROBE.SimID." + to_string(Settings::SimID) + ".MaxSepClustSzDistr."
			+ whichtarget + "." + againstwhich + ".TskID." + to_string(tid) + ".RunID." + to_string(runid) + ".csv";

	ofstream sslog;
	sslog.open( fn.c_str() );
	if ( sslog.is_open() == false ) {
		string mess = "Unable to open file " + fn;
		stopping(mess);
		return;
	}

	sslog.precision(18);
	sslog << "ClusterID;NumberOfMembers;BaryCenterX;BaryCenterY;BaryCenterZ\n";
	sslog << "1;1;nm;nm;nm\n";
	sslog << "ClusterID;NumberOfMembers;BaryCenterX;BaryCenterY;BaryCenterZ\n";
	size_t excluded = 0;
	for( auto it = precipitates.begin(); it != precipitates.end(); it++ ) {
		if ( it->has_boundary_contact() == false ) {
			p3d bary = it->center_of_gravity();
			sslog << it->id << ";" << it->n << ";" << bary.x << ";" << bary.y << ";" << bary.z << "\n";
		}
		else {
			excluded++;
		}
	}
	sslog.flush();
	sslog.close();

	cout << "A total of " << excluded << " out of " << precipitates.size() << " cluster were expelled from the size distribution because of boundary contact!" << endl;
}

void clusterHdl::report_size_distr_vtk( const string whichtarget, const string againstwhich, const unsigned int tid, const unsigned int runid )
{
	//MK::write VTK file showing barycenter of all synthesized cluster in reconstructed space
	//includes particles outside actual tip

	string vtk_io_fn = "PARAPROBE.SimID." + to_string(Settings::SimID) + ".MaxSepClustSzDistr."
				+ whichtarget + "." + againstwhich + ".TskID." + to_string(tid) + ".RunID." + to_string(runid) + ".vtk";

	ofstream vtk;
	vtk.open(vtk_io_fn.c_str(), ofstream::out | ofstream::trunc);
	if (vtk.is_open() == true) {
		//header
		vtk << "# vtk DataFile Version 2.0\n";
		vtk << "PARAPROBE HPDBScan clustering " << whichtarget << " " << againstwhich << " " << tid << " " << " runid " << runid << "\n";
		vtk << "ASCII\n";
		vtk << "DATASET POLYDATA\n";
		vtk << "POINTS " << precipitates.size() << " double\n";
		for( auto it = precipitates.begin(); it != precipitates.end(); it++ ) {
			p3d bary = it->center_of_gravity();
			vtk << bary.x << " " << bary.y << " " << bary.z << "\n";
		}
		vtk << "\n";
		vtk << "VERTICES " << precipitates.size() << " " << 2*precipitates.size() << "\n";
		for( size_t i = 0; i < precipitates.size(); ++i ) {
			vtk << 1 << " " << i << "\n";
		}
		vtk << "\n";
		vtk << "POINT_DATA " << precipitates.size() << "\n";
		vtk << "FIELD FieldData 1\n";
		vtk << "NumberOfMembers 1 " << precipitates.size() << " float\n";
		for( auto it = precipitates.begin(); it != precipitates.end(); it++ ) {
			vtk << it->n << "\n";
		}
		vtk.flush();
		vtk.close();
		cout << "VTK file of clustered particle positions written to file" << endl;
	}
	else {
		cout << "VTK file of clustered particle positions was not written" << endl;
	}
}
