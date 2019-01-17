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


#include "PARAPROBE_AABBTree.h"

/*
  Copyright (c) 2009 Erin Catto http://www.box2d.org
  Copyright (c) 2016-2017 Lester Hedges <lester.hedges+aabbcc@gmail.com>

  This software is provided 'as-is', without any express or implied
  warranty. In no event will the authors be held liable for any damages
  arising from the use of this software.

  Permission is granted to anyone to use this software for any purpose,
  including commercial applications, and to alter it and redistribute it
  freely, subject to the following restrictions:

  1. The origin of this software must not be misrepresented; you must not
     claim that you wrote the original software. If you use this software
     in a product, an acknowledgment in the product documentation would be
     appreciated but is not required.

  2. Altered source versions must be plainly marked as such, and must not be
     misrepresented as being the original software.

  3. This notice may not be removed or altered from any source distribution.

  This code was adapted from parts of the Box2D Physics Engine,
  http://www.box2d.org
*/

//MK::modifications considering that the apt_real world APT tip is 3d and does not show periodicity

AABB::AABB(const trpl lowerBound_, const trpl upperBound_)
{
    lowerBound = lowerBound_;
	upperBound = upperBound_;

	surfaceArea = computeSurfaceArea();
    centre = computeCentre();
}

apt_real AABB::computeSurfaceArea() const
{
	// Calculate the surface area of the 3D AABB.
	apt_real wx = upperBound.a - lowerBound.a;
	apt_real wy = upperBound.b - lowerBound.b;
	apt_real wz = upperBound.c - lowerBound.c;
	return 2.f * (wx*wy + wx*wz + wy*wz);
}

apt_real AABB::getSurfaceArea() const
{
    return surfaceArea;
}

void AABB::merge(const AABB& aabb1, const AABB& aabb2)
{
	//assert not required because aabb1 and aabb2 both 3d
	lowerBound.a = min( aabb1.lowerBound.a, aabb2.lowerBound.a );
	lowerBound.b = min( aabb1.lowerBound.b, aabb2.lowerBound.b );
	lowerBound.c = min( aabb1.lowerBound.c, aabb2.lowerBound.c );

	upperBound.a = max( aabb1.upperBound.a, aabb2.upperBound.a );
	upperBound.b = max( aabb1.upperBound.b, aabb2.upperBound.b );
	upperBound.c = max( aabb1.upperBound.c, aabb2.upperBound.c );

    surfaceArea = computeSurfaceArea();
    centre = computeCentre();
}

bool AABB::contains(const AABB& aabb) const
{
	//assert not required because aabb and myself both 3d
	if (lowerBound.a < aabb.lowerBound.a) return false;
	if (upperBound.a > aabb.upperBound.a) return false;
 	if (lowerBound.b < aabb.lowerBound.b) return false;
	if (upperBound.b > aabb.upperBound.b) return false; 
	if (lowerBound.c < aabb.lowerBound.c) return false;
	if (upperBound.c > aabb.upperBound.c) return false;
	return true;
}

bool AABB::overlaps(const AABB& aabb) const
{
	//assert not required because aabb and myself both 3d
	//MK::only 3d case to evaluate
	return !(
		aabb.upperBound.a < lowerBound.a 
		|| aabb.lowerBound.a > upperBound.a
		|| aabb.upperBound.b < lowerBound.b
		|| aabb.lowerBound.b > upperBound.b
		|| aabb.upperBound.c < lowerBound.c
		|| aabb.lowerBound.c > upperBound.c
	);
}

trpl AABB::computeCentre()
{
    return trpl(0.5*(lowerBound.a+upperBound.a), 0.5*(lowerBound.b+upperBound.b), 0.5*(lowerBound.c+upperBound.c) );
}



bool Node::isLeaf() const
{
    return (left == NULL_NODE);
}

Tree::Tree(unsigned int nParticles)
{
    // Initialise the tree.
    root = NULL_NODE;
    nodeCount = 0;
    nodeCapacity = nParticles;
    nodes.resize(nodeCapacity);

    // Build a linked list for the list of free nodes.
    for (unsigned int i=0;i<nodeCapacity-1;i++)
    {
        nodes[i].next = i + 1;
        nodes[i].height = -1;
    }
    nodes[nodeCapacity-1].next = NULL_NODE;
    nodes[nodeCapacity-1].height = -1;

    // Assign the index of the first free node.
    freeList = 0;
}


unsigned int Tree::allocateNode()
{
	// Exand the node pool as needed.
	if (freeList == NULL_NODE) {
		assert(nodeCount == nodeCapacity);
		//if ( nodeCount == nodeCapacity ) {
			// The free list is empty. Rebuild a bigger pool.

			//##MK::original quick-and-dirty unsigned int * 2 may overflow!
			unsigned long safeCapacityIncr = static_cast<unsigned long>(nodeCapacity);
			safeCapacityIncr *= 2;
			if ( safeCapacityIncr < static_cast<unsigned long>(numeric_limits<unsigned int>::max() - 2) ) { //will now not overflow on nodeCapacity
				nodeCapacity *= 2;
				nodes.resize(nodeCapacity);
			}
			else {
				//would at least touch flag ID or even overflow therefor reject
				return TREE_ERROR;
			}

			// Build a linked list for the list of free nodes.
			for (unsigned int i=nodeCount; i<nodeCapacity-1; i++) {
				nodes[i].next = i + 1;
				nodes[i].height = -1;
			}
			nodes[nodeCapacity-1].next = NULL_NODE;
			nodes[nodeCapacity-1].height = -1;

			// Assign the index of the first free node.
			freeList = nodeCount;
		//}
		//else { //on error return NULL_NODE;
		//	return TREE_ERROR;
		//}
	}

	// Peel a node off the free list.
	unsigned int node = freeList;
	freeList = nodes[node].next;
	nodes[node].parent = NULL_NODE;
	nodes[node].left = NULL_NODE;
	nodes[node].right = NULL_NODE;
	nodes[node].height = 0;
	nodeCount++;

	return node;
}


void Tree::freeNode(unsigned int node)
{
    //assert( 0 <= node && node < nodeCapacity );
	//assert( 0 < nodeCount ); 
	//does not make much sense because by type node is positive or zero anyway and nodeCount <= nodeCapacity < numeric_limits<unsigned int>::max()-1
	
	assert( node < nodeCount );

	//if ( node < nodeCount ) {
		nodes[node].next = freeList;
		nodes[node].height = -1;
		freeList = node;
		nodeCount--;
		//return 0;
	//}
	//else {
	//	return TREE_ERROR;
	//}
}


void Tree::insertParticle(unsigned int particle, const trpl position, apt_real radius)
{
    // Allocate a new node for the particle.
    unsigned int node = allocateNode();

	nodes[node].aabb.lowerBound.a = position.a - radius;
	nodes[node].aabb.upperBound.a = position.a + radius;
	nodes[node].aabb.lowerBound.b = position.b - radius;
	nodes[node].aabb.upperBound.b = position.b + radius;
	nodes[node].aabb.lowerBound.c = position.c - radius;
	nodes[node].aabb.upperBound.c = position.c + radius;

    nodes[node].aabb.surfaceArea = nodes[node].aabb.computeSurfaceArea();
    nodes[node].aabb.centre = nodes[node].aabb.computeCentre();

    // Zero the height.
    nodes[node].height = 0;

    // Insert a new leaf into the tree.
	if ( insertLeaf(node) == TREE_ERROR ) {
		cout << "ERROR::Tree leaf insertion error!" << endl;
	}

    // Add the new particle to the map.
    particleMap.insert(map<unsigned int, unsigned int>::value_type(particle, node));

    // Store the particle index.
    nodes[node].particle = particle;
}

void Tree::insertParticle(unsigned int particle, const trpl lowerBound, const trpl upperBound)
{
    // Allocate a new node for the particle.
    unsigned int node = allocateNode();

    // Compute the AABB limits
	nodes[node].aabb.lowerBound = lowerBound;
	nodes[node].aabb.upperBound = upperBound;
	    
	//MK::no fattening needed as measured locations in tip are static
    nodes[node].aabb.surfaceArea = nodes[node].aabb.computeSurfaceArea();
    nodes[node].aabb.centre = nodes[node].aabb.computeCentre();

    // Zero the height.
    nodes[node].height = 0;

    // Insert a new leaf into the tree.
	if ( insertLeaf(node) == TREE_ERROR ) {
		cout << "ERROR::Tree leaf insertion error!" << endl;
	}

    // Add the new particle to the map.
    particleMap.insert(map<unsigned int, unsigned int>::value_type(particle, node));

    // Store the particle index.
    nodes[node].particle = particle;
}

void Tree::removeParticle(unsigned int particle)
{
    // Map iterator.
    map<unsigned int, unsigned int>::iterator it;

    // Find the particle.
    it = particleMap.find(particle);

    // Extract the node index.
    unsigned int node = it->second;

    // Erase the particle from the map.
    particleMap.erase(it);

	//assert( 0 <= node && node < nodeCapacity ); //does not make much sense because node is unsigned int anyway...
	assert( node < nodeCapacity );
	assert( nodes[node].isLeaf() );

	//unsigned int status = 0;
    //if ( node < nodeCapacity && nodes[node].isLeaf() == true ) { //assert(0 <= node && node < nodeCapacity); //assert(nodes[node].isLeaf());
		removeLeaf(node);
		//status = freeNode(node);
		freeNode(node);
	//}
	//else { //error
	//	status = numeric_limits<unsigned int>::max();
	//}	
}


vector<unsigned int> Tree::query(unsigned int particle)
{
    if ( particleMap.count(particle) > 0 ) {
		// Test overlap of particle AABB against all other particles.
		return query(particle, nodes[particleMap.find(particle)->second].aabb);
	}
	else {
		vector<unsigned int> nothing;
		return nothing;
	}
}


vector<unsigned int> Tree::query(unsigned int particle, const AABB& aabb)
{
    vector<unsigned int> stack;
    stack.reserve(256);
    stack.push_back(root);

    vector<unsigned int> particles;

    while (stack.size() > 0)
    {
        unsigned int node = stack.back();
        stack.pop_back();

        // Copy the AABB.
        AABB nodeAABB = nodes[node].aabb;

        if (node == NULL_NODE) continue;

        // Test for overlap between the AABBs.
        if (aabb.overlaps(nodeAABB))
        {
            // Check that we're at a leaf node.
            if (nodes[node].isLeaf())
            {
                // Can't interact with itself.
                if (nodes[node].particle != particle)
                    particles.push_back(nodes[node].particle);
            }
            else
            {
                stack.push_back(nodes[node].left);
                stack.push_back(nodes[node].right);
            }
        }
    }

    return particles;
}

vector<unsigned int> Tree::query(const AABB& aabb)
{
    // Test overlap of AABB against all particles.
    return query(numeric_limits<unsigned int>::max(), aabb);
}


const AABB& Tree::getAABB(unsigned int particle)
{
    return nodes[particleMap[particle]].aabb;
}


unsigned int Tree::insertLeaf(unsigned int leaf)
{
    if (root == NULL_NODE)
    {
        root = leaf;
        nodes[root].parent = NULL_NODE;
        return 0;
    }

    // Find the best sibling for the node.

    AABB leafAABB = nodes[leaf].aabb;
    unsigned int index = root;

    while (!nodes[index].isLeaf())
    {
        // Extract the children of the node.
        unsigned int left  = nodes[index].left;
        unsigned int right = nodes[index].right;

        apt_real surfaceArea = nodes[index].aabb.getSurfaceArea();

        AABB combinedAABB;
        combinedAABB.merge(nodes[index].aabb, leafAABB);
        apt_real combinedSurfaceArea = combinedAABB.getSurfaceArea();

        // Cost of creating a new parent for this node and the new leaf.
        apt_real cost = 2.0 * combinedSurfaceArea;

        // Minimum cost of pushing the leaf further down the tree.
        apt_real inheritanceCost = 2.0 * (combinedSurfaceArea - surfaceArea);

        // Cost of descending to the left.
        apt_real costLeft;
        if (nodes[left].isLeaf())
        {
            AABB aabb;
            aabb.merge(leafAABB, nodes[left].aabb);
            costLeft = aabb.getSurfaceArea() + inheritanceCost;
        }
        else
        {
            AABB aabb;
            aabb.merge(leafAABB, nodes[left].aabb);
            apt_real oldArea = nodes[left].aabb.getSurfaceArea();
            apt_real newArea = aabb.getSurfaceArea();
            costLeft = (newArea - oldArea) + inheritanceCost;
        }

        // Cost of descending to the right.
        apt_real costRight;
        if (nodes[right].isLeaf())
        {
            AABB aabb;
            aabb.merge(leafAABB, nodes[right].aabb);
            costRight = aabb.getSurfaceArea() + inheritanceCost;
        }
        else
        {
            AABB aabb;
            aabb.merge(leafAABB, nodes[right].aabb);
            apt_real oldArea = nodes[right].aabb.getSurfaceArea();
            apt_real newArea = aabb.getSurfaceArea();
            costRight = (newArea - oldArea) + inheritanceCost;
        }

        // Descend according to the minimum cost.
        if ((cost < costLeft) && (cost < costRight)) break;

        // Descend.
        if (costLeft < costRight) index = left;
        else                      index = right;
    }

    unsigned int sibling = index;

    // Create a new parent.
    unsigned int oldParent = nodes[sibling].parent;
    unsigned int newParent = allocateNode();
    nodes[newParent].parent = oldParent;
    nodes[newParent].aabb.merge(leafAABB, nodes[sibling].aabb);
    nodes[newParent].height = nodes[sibling].height + 1;

    // The sibling was not the root.
    if (oldParent != NULL_NODE)
    {
        if (nodes[oldParent].left == sibling) nodes[oldParent].left = newParent;
        else                                  nodes[oldParent].right = newParent;

        nodes[newParent].left = sibling;
        nodes[newParent].right = leaf;
        nodes[sibling].parent = newParent;
        nodes[leaf].parent = newParent;
    }
    // The sibling was the root.
    else
    {
        nodes[newParent].left = sibling;
        nodes[newParent].right = leaf;
        nodes[sibling].parent = newParent;
        nodes[leaf].parent = newParent;
        root = newParent;
    }

    // Walk back up the tree fixing heights and AABBs.
    index = nodes[leaf].parent;
    while (index != NULL_NODE)
    {
        index = balance(index);

        unsigned int left = nodes[index].left;
        unsigned int right = nodes[index].right;

		assert( left != NULL_NODE );
		assert( right != NULL_NODE );

        //if ( left != NULL_NODE && right != NULL_NODE ) {

			nodes[index].height = 1 + max(nodes[left].height, nodes[right].height);
			nodes[index].aabb.merge(nodes[left].aabb, nodes[right].aabb);

			index = nodes[index].parent;
		//}
		//else { return TREE_ERROR; }
    }

	return 0;
}

void Tree::removeLeaf(unsigned int leaf)
{
    if (leaf == root)
    {
        root = NULL_NODE;
        return;
    }

    unsigned int parent = nodes[leaf].parent;
    unsigned int grandParent = nodes[parent].parent;
    unsigned int sibling;

    if (nodes[parent].left == leaf) sibling = nodes[parent].right;
    else                            sibling = nodes[parent].left;

    // Destroy the parent and connect the sibling to the grandparent.
    if (grandParent != NULL_NODE)
    {
        if (nodes[grandParent].left == parent) nodes[grandParent].left = sibling;
        else                                   nodes[grandParent].right = sibling;

        nodes[sibling].parent = grandParent;
        //unsigned int status = freeNode(parent);
		freeNode(parent);

        // Adjust ancestor bounds.
        unsigned int index = grandParent;
        while (index != NULL_NODE)
        {
            index = balance(index);

            unsigned int left = nodes[index].left;
            unsigned int right = nodes[index].right;

            nodes[index].aabb.merge(nodes[left].aabb, nodes[right].aabb);
            nodes[index].height = 1 + max(nodes[left].height, nodes[right].height);

            index = nodes[index].parent;
        }
    }
    else
    {
        root = sibling;
        nodes[sibling].parent = NULL_NODE;
        //unsigned int status = freeNode(parent);
		freeNode(parent);
    }
}

unsigned int Tree::balance(unsigned int node)
{
    assert(node != NULL_NODE);

    if (nodes[node].isLeaf() || (nodes[node].height < 2))
        return node;

    unsigned int left = nodes[node].left;
    unsigned int right = nodes[node].right;

    //assert((0 <= left) && (left < nodeCapacity));
    //assert((0 <= right) && (right < nodeCapacity));
	assert( left < nodeCapacity );
	assert( right < nodeCapacity );

    int currentBalance = nodes[right].height - nodes[left].height;

    // Rotate right branch up.
    if (currentBalance > 1)
    {
        unsigned int rightLeft = nodes[right].left;
        unsigned int rightRight = nodes[right].right;

		//these asserts are unnecessarily inefficient as right** are unsigned ints thus 0>= implicitly!
        //assert((0 <= rightLeft) && (rightLeft < nodeCapacity));
        //assert((0 <= rightRight) && (rightRight < nodeCapacity));
		assert( rightLeft < nodeCapacity );
		assert( rightRight < nodeCapacity );

        // Swap node and its right-hand child.
        nodes[right].left = node;
        nodes[right].parent = nodes[node].parent;
        nodes[node].parent = right;

        // The node's old parent should now point to its right-hand child.
        if (nodes[right].parent != NULL_NODE)
        {
            if (nodes[nodes[right].parent].left == node) 
				nodes[nodes[right].parent].left = right;
            else
            {
                assert(nodes[nodes[right].parent].right == node);
                nodes[nodes[right].parent].right = right;
            }
        }
        else root = right;

        // Rotate.
        if (nodes[rightLeft].height > nodes[rightRight].height)
        {
            nodes[right].right = rightLeft;
            nodes[node].right = rightRight;
            nodes[rightRight].parent = node;
            nodes[node].aabb.merge(nodes[left].aabb, nodes[rightRight].aabb);
            nodes[right].aabb.merge(nodes[node].aabb, nodes[rightLeft].aabb);

            nodes[node].height = 1 + max(nodes[left].height, nodes[rightRight].height);
            nodes[right].height = 1 + max(nodes[node].height, nodes[rightLeft].height);
        }
        else
        {
            nodes[right].right = rightRight;
            nodes[node].right = rightLeft;
            nodes[rightLeft].parent = node;
            nodes[node].aabb.merge(nodes[left].aabb, nodes[rightLeft].aabb);
            nodes[right].aabb.merge(nodes[node].aabb, nodes[rightRight].aabb);

            nodes[node].height = 1 + max(nodes[left].height, nodes[rightLeft].height);
            nodes[right].height = 1 + max(nodes[node].height, nodes[rightRight].height);
        }

        return right;
    }

    // Rotate left branch up.
    if (currentBalance < -1)
    {
        unsigned int leftLeft = nodes[left].left;
        unsigned int leftRight = nodes[left].right;

        //assert((0 <= leftLeft) && (leftLeft < nodeCapacity));
        //assert((0 <= leftRight) && (leftRight < nodeCapacity));
		assert(leftLeft < nodeCapacity);
		assert(leftRight < nodeCapacity);

        // Swap node and its left-hand child.
        nodes[left].left = node;
        nodes[left].parent = nodes[node].parent;
        nodes[node].parent = left;

        // The node's old parent should now point to its left-hand child.
        if (nodes[left].parent != NULL_NODE)
        {
            if (nodes[nodes[left].parent].left == node) 
				nodes[nodes[left].parent].left = left;
            else
            {
                assert(nodes[nodes[left].parent].right == node);
                nodes[nodes[left].parent].right = left;
            }
        }
        else root = left;

        // Rotate.
        if (nodes[leftLeft].height > nodes[leftRight].height)
        {
            nodes[left].right = leftLeft;
            nodes[node].left = leftRight;
            nodes[leftRight].parent = node;
            nodes[node].aabb.merge(nodes[right].aabb, nodes[leftRight].aabb);
            nodes[left].aabb.merge(nodes[node].aabb, nodes[leftLeft].aabb);

            nodes[node].height = 1 + max(nodes[right].height, nodes[leftRight].height);
            nodes[left].height = 1 + max(nodes[node].height, nodes[leftLeft].height);
        }
        else
        {
            nodes[left].right = leftRight;
            nodes[node].left = leftLeft;
            nodes[leftLeft].parent = node;
            nodes[node].aabb.merge(nodes[right].aabb, nodes[leftLeft].aabb);
            nodes[left].aabb.merge(nodes[node].aabb, nodes[leftRight].aabb);

            nodes[node].height = 1 + max(nodes[right].height, nodes[leftLeft].height);
            nodes[left].height = 1 + max(nodes[node].height, nodes[leftRight].height);
        }

        return left;
    }

    return node;
}

unsigned int Tree::computeHeight() const
{
    return computeHeight(root);
}

unsigned int Tree::computeHeight(unsigned int node) const
{
    //assert((0 <= node) && (node < nodeCapacity));
	assert(node < nodeCapacity);

    if (nodes[node].isLeaf()) return 0;

    unsigned int height1 = computeHeight(nodes[node].left);
    unsigned int height2 = computeHeight(nodes[node].right);

    return 1 + max(height1, height2);
}

unsigned int Tree::getHeight() const
{
    if (root == NULL_NODE) 
		return 0;
    return nodes[root].height;
}

unsigned int Tree::getNodeCount() const
{
    return nodeCount;
}

unsigned int Tree::computeMaximumBalance() const
{
    unsigned int maxBalance = 0;
    for (unsigned int i=0; i<nodeCapacity; i++)
    {
        if (nodes[i].height <= 1)
            continue;

        assert(nodes[i].isLeaf() == false);

        unsigned int balance = abs(nodes[nodes[i].left].height - nodes[nodes[i].right].height);
        maxBalance = max(maxBalance, balance);
    }

    return maxBalance;
}

apt_real Tree::computeSurfaceAreaRatio() const
{
    if (root == NULL_NODE)
    	return 0.0;

    apt_real rootArea = nodes[root].aabb.computeSurfaceArea();
    apt_real totalArea = 0.0;

    for (unsigned int i=0; i<nodeCapacity;i++)
    {
        if (nodes[i].height < 0) continue;

        totalArea += nodes[i].aabb.computeSurfaceArea();
    }

    return totalArea / rootArea;
}


void Tree::validate() const
{
	#ifndef NDEBUG
		validateStructure(root);
		validateMetrics(root);

		unsigned int freeCount = 0;
		unsigned int freeIndex = freeList;
		while (freeIndex != NULL_NODE) {
			//assert((0 <= freeIndex) && (freeIndex < nodeCapacity));
			assert( freeIndex < nodeCapacity);

			freeIndex = nodes[freeIndex].next;
			freeCount++;
		}
		assert(getHeight() == computeHeight());
		assert((nodeCount + freeCount) == nodeCapacity);
	#endif
}


void Tree::report_tree( const string csv_fn )
{
	ofstream csvlog;
	csvlog.open( csv_fn.c_str() );
	if ( csvlog.is_open() == false ) {
		string mess = "Unable to open file " + csv_fn;
		stopping(mess);
		return;
	}

	csvlog.precision(18);
	csvlog << "NodeID;XMin;YMin;ZMin;XMax;YMax;ZMax;IsLeaf\n";
	csvlog << "1;nm;nm;nm;nm;nm;nm;bool\n";
	csvlog << "NodeID;XMin;YMin;ZMin;XMax;YMax;ZMax;IsLeaf\n";

	size_t nodeidx = 0;
	for( auto it = nodes.begin(); it != nodes.end(); ++it, ++nodeidx ) {
		csvlog << nodeidx << ";";
		csvlog << it->aabb.lowerBound.a << ";" << it->aabb.lowerBound.b << ";" << it->aabb.lowerBound.c << ";";
		csvlog << it->aabb.upperBound.a << ";" << it->aabb.upperBound.b << ";" << it->aabb.upperBound.c << ";";
		csvlog << it->isLeaf() << "\n";
	}
	csvlog.flush();
	csvlog.close();
}


void Tree::validateStructure(unsigned int node) const
{
	if (node == NULL_NODE) 
		return;
	
	if (node == root) 
		assert(nodes[node].parent == NULL_NODE);

	unsigned int left = nodes[node].left;
	unsigned int right = nodes[node].right;

	if (nodes[node].isLeaf() ) {
		assert(left == NULL_NODE);
		assert(right == NULL_NODE);
		assert(nodes[node].height == 0);
		return;
	}
	//assert((0 <= left) && (left < nodeCapacity));
	//assert((0 <= right) && (right < nodeCapacity));
	assert( left < nodeCapacity);
	assert( right < nodeCapacity);	
	assert(nodes[left].parent == node);
	assert(nodes[right].parent == node);

	validateStructure(left);
	validateStructure(right);
}


void Tree::validateMetrics(unsigned int node) const
{
	if (node == NULL_NODE) 
		return;

	unsigned int left = nodes[node].left;
	unsigned int right = nodes[node].right;

	if (nodes[node].isLeaf() ) {
		assert(left == NULL_NODE);
		assert(right == NULL_NODE);
		assert(nodes[node].height == 0);
		return;
	}
	//assert((0 <= left) && (left < nodeCapacity));
	//assert((0 <= right) && (right < nodeCapacity));
	assert( left < nodeCapacity);
	assert( right < nodeCapacity);

	int height1 = nodes[left].height;
	int height2 = nodes[right].height;
	int height = 1 + max(height1, height2);
	assert(nodes[node].height == height);

	AABB aabb;
	aabb.merge(nodes[left].aabb, nodes[right].aabb);

	assert( (aabb.lowerBound.a - nodes[node].aabb.lowerBound.a) < ASSERTATION_ACCURACY );
	assert( (aabb.upperBound.a - nodes[node].aabb.upperBound.a) < ASSERTATION_ACCURACY );
	assert( (aabb.lowerBound.b - nodes[node].aabb.lowerBound.b) < ASSERTATION_ACCURACY );
	assert( (aabb.upperBound.b - nodes[node].aabb.upperBound.b) < ASSERTATION_ACCURACY );
	assert( (aabb.lowerBound.c - nodes[node].aabb.lowerBound.c) < ASSERTATION_ACCURACY );
	assert( (aabb.upperBound.c - nodes[node].aabb.upperBound.c) < ASSERTATION_ACCURACY );

	validateMetrics(left);
	validateMetrics(right);
}
