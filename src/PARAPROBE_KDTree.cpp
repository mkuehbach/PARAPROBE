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

#include "PARAPROBE_KDTree.h"


void kd_tree::build(const vector<p3d> & points, vector<size_t> & permutate )
{
	//const size_t page_size = 4096;
	const size_t leaf_size = 256; //page_size / sizeof(p3d);

//cout << "sizeof(p3dm1) = " << sizeof(p3d) << endl;
//cout << "Leaf size = " << leaf_size << endl;

	const size_t stack_size = 32;
	const size_t count = points.size();
	if (count == 0)
		return;

	for (size_t i = 0; i < count; ++i)
		permutate.push_back(i);

	nodes.push_back( node() );
	build_task tasks[stack_size] = {0, count-1, 0, 0}; //accessing on i0, i1, [first, last], node id, split dim
	int current_task = 0;
	p3d min(RMAX, RMAX, RMAX);
	p3d max(RMIN, RMIN, RMIN);

	do
	{
		build_task task = tasks[current_task];
		node &n = nodes[task.node];
		if ( (task.last - task.first) < leaf_size ) {
			//leaf_threshold, a window of the permutate field that is sorted [task.first <= task.last]
			//n.split MUST NOT BE SET! as it indicates that this node specifies now a leaf, and therefore i0 and i1 reference positions on permutation field rather on this->nodes !
			//##MK::DEBUG likewise we leave n.dimdebug -1
			//likewise n.splitpos needs to be numeric_limits<apt_xyz>::max()
			//assert( n.dimdebug == -1);
			assert( n.splitpos == numeric_limits<apt_xyz>::max() );


			//i.e. a sorted contiguous region on the permutation index field whose indices specify the points in this leaf
			//in case that a node is a leaf ONLY, i0 and i1 specify task.first and task.last as markers where we read indices from the permutation field index inclusively
			n.i0 = task.first; //MK::now transforming ID index i0 from indexing nodes to points
			n.i1 = task.last;
			//n.i0 and n.i1 can be the same namely if only single in the leaf
			for(size_t k = task.first; k <= task.last; ++k) {
				expand(min, max, points[permutate[k]]);
			}

			assert(is_leaf(n));
			--current_task;
			continue;
		}

		//no leaf yet so split
		const size_t k = (task.first + task.last) / 2;
//cout << "k/first/last = " << k << "\t\t" << task.first << "\t\t" << task.last << endl;
		//n.dimdebug = task.dim;

		if ( task.dim == PARAPROBE_XAXIS) {
			nth_element( permutate.begin() + task.first, permutate.begin() + k, permutate.begin() + task.last + 1, lower_x(points) );
			n.splitpos = points[permutate[k]].x;
		}
		else if ( task.dim == PARAPROBE_YAXIS ) {
			nth_element( permutate.begin() + task.first, permutate.begin() + k, permutate.begin() + task.last + 1, lower_y(points) );
			n.splitpos = points[permutate[k]].y;
		}
		else { //task.dim == PARAPROBE_ZAXIS
			nth_element( permutate.begin() + task.first, permutate.begin() + k, permutate.begin() + task.last + 1, lower_z(points) );
			n.splitpos = points[permutate[k]].z;
		}

		const size_t i0 = nodes.size(); //MK::i0 and i1 here indexing nodes, not points!
		const size_t i1 = i0 + 1;

		n.split = k;

		n.i0 = i0;
		n.i1 = i1;

		size_t next_dim = ((task.dim+1) > PARAPROBE_ZAXIS) ? PARAPROBE_XAXIS : (task.dim + 1); //task.dim++; //cyclic X,Y,Z, X, Y, ...

		const build_task taskleft = {task.first, k, i0, next_dim};
		const build_task taskright = {k + 1, task.last, i1, next_dim};
		tasks[current_task] = taskleft;
		nodes.push_back( node() );
		++current_task;

		assert( static_cast<size_t>(current_task) < stack_size);

		tasks[current_task] = taskright;
		//do not add ++current_task because we replace the existent on the stack on which we still work "so two new for an old"
		nodes.push_back( node() );
	} while ( current_task != -1 );


	//#pragma omp critical
	//{
	//	cout << "KDTree resident memory at least " << get_treememory_consumption() << " B" << endl;
	//	cout << "Nodes in total " << nodes.size() << endl;
	//}


	this->min = min;
	this->max = max;
}

void kd_tree::pack_p3dm1_d( const vector<size_t> & permutate, const vector<p3dm1> & apt1, const vector<apt_xyz> & dist1, vector<p3dm1> & apt2, vector<apt_xyz> & dist2 )
{
	//permutate array holds indices of points from apt1 that are often randomly shuffled
	//to avoid future cache inefficiency when querying we fuse and arrange linearly in memory according to tree node ownership
	//therewith SIMD on contiguous memory blocks can be applied more efficiently at the leaf level

cout << "KDTree build complete reorganizing input array" << endl;
	const size_t count = apt1.size();
	assert( count == permutate.size() );
	assert( count == dist1.size() );
	for(size_t i = 0; i < count; ++i) { //MK::also random memory access
		p3dm1 thispoint = apt1[permutate[i]];
		apt_xyz thisdistance = dist1[permutate[i]];
		apt2.push_back( thispoint );
		dist2.push_back( thisdistance );
	}
}

void kd_tree::pack_p3dm1( const vector<size_t> & permutate, const vector<p3d> & aptpos1, const vector<unsigned int> & aptlabel1, vector<p3dm1> & apt2 )
{
	//permutate array holds indices of points from aptpos1 and corresponding ion range labels aptlabel1 that are often randomly shuffled
	//to avoid future cache inefficiency when querying we fuse the individual pieces of information in aptpos1 and aptlabel1 and
	//arrange linear in memory according to tree node ownership
	//therewith SIMD on contiguous memory blocks can be applied more efficiently at the leaf level

cout << "KDTree build complete reorganizing input array" << endl;
	const size_t count = aptpos1.size();
	assert( count == aptlabel1.size() );
	assert( count == permutate.size() );
	for(size_t i = 0; i < count; ++i) { //MK::also random memory access
		p3d thispoint = aptpos1[permutate[i]];
		unsigned int thislabel = aptlabel1[permutate[i]];
		apt2.push_back( p3dm1(thispoint.x, thispoint.y, thispoint.z, thislabel) );
	}
}


void kd_tree::pack_p3dm1_dist( vector<size_t> const & permutate, vector<p3dm1> const & in1, vector<apt_xyz> const & in2, vector<p3dm1> & out1, vector<apt_xyz> & out2 )
{
	//permutate array holds indices of points from aptpos1 and corresponding ion range labels aptlabel1 that are often randomly shuffled
	//to avoid future cache inefficiency when querying we fuse the individual pieces of information in aptpos1 and aptlabel1 and
	//arrange linear in memory according to tree node ownership
	//therewith SIMD on contiguous memory blocks can be applied more efficiently at the leaf level

cout << "KDTree build complete reorganizing input array" << endl;
	const size_t count = in1.size();
	assert( count == in2.size() );
	assert( count == permutate.size() );
	for(size_t i = 0; i < count; ++i) { //MK::also random memory access
		p3dm1 thismarkedpoint = in1[permutate[i]];
		apt_xyz thisdistance = in2[permutate[i]]; //##MK::conflict of interest two cache misses per entry, better build KDTree as first step

		out1.push_back( thismarkedpoint );
		out2.push_back( thisdistance );
	}
}



p3dm1 kd_tree::nearest(const size_t idx, const vector<p3dm1> & sortedpoints, apt_xyz epsball ) const
{
	const size_t stack_size = 32;
	const cuboid aabb(min, max);
	p3dm1 target = sortedpoints[idx];
	apt_xyz currentbestd = epsball; //MK::squared value to avoid sqrt computations

	traverse_task tasks[stack_size] = {aabb, 0, 0, aabb.outside_proximity(target)}; //first box, first node, dim, distance to first
	int current_task = 0;

	p3dm1 best = p3dm1();
	do
	{
		const traverse_task &task = tasks[current_task];
		if (task.d > currentbestd) {
			--current_task;
			continue;
		}

		const node &n = nodes[task.node];
		if ( is_leaf(n) == true ) {
			//if n is a node n.i0 and n.i1 specify contiguous range of points on SORTED!! sortedpoints
			for(size_t i = n.i0; i <= n.i1; ++i) {
				//##MK::do SIMD
				const p3dm1 query = sortedpoints[i];
				const apt_xyz d = euclidean_sqrd(target, query);
				if ( d > currentbestd ) {
					continue;
				} else {
					if ( idx != i ) {
						currentbestd = d;
						best = query;
					}
				}
			}
			--current_task;
			continue;
		}

		traverse_task near;
		traverse_task far;
		//assert( n.dimdebug == task.dim );
		if ( task.dim == PARAPROBE_XAXIS ) {
			const apt_xyz splitposx = n.splitpos; //MK::using node-local piece of information rather than probing random memory position
			const cuboid left = cuboid( task.bx.min, p3d(splitposx, task.bx.max.y, task.bx.max.z) ); //#####MK::z coordinate of point...
			const cuboid right = cuboid( p3d(splitposx, task.bx.min.y, task.bx.min.z), task.bx.max );
			if ( target.x <= splitposx ) {
				near.bx = left;
				near.node = n.i0;
				far.bx = right;
				far.node = n.i1;
			}
			else {
				near.bx = right;
				near.node = n.i1;
				far.bx = left;
				far.node = n.i0;
			}
		}
		else if ( task.dim == PARAPROBE_YAXIS ) {
			const apt_xyz splitposy = n.splitpos;
			const cuboid lower = cuboid( task.bx.min, p3d(task.bx.max.x, splitposy, task.bx.max.z) );
			const cuboid upper = cuboid( p3d(task.bx.min.x, splitposy, task.bx.min.z), task.bx.max );
			if ( target.y <= splitposy ) {
				near.bx = lower;
				near.node = n.i0;
				far.bx = upper;
				far.node = n.i1;
			}
			else {
				near.bx = upper;
				near.node = n.i1;
				far.bx = lower;
				far.node = n.i0;
			}
		}
		else {
			const apt_xyz splitposz = n.splitpos;
			const cuboid back = cuboid( task.bx.min, p3d(task.bx.max.x, task.bx.max.y, splitposz) );
			const cuboid front = cuboid( p3d(task.bx.min.x, task.bx.min.y, splitposz), task.bx.max );
			if ( target.z <= splitposz ) {
				near.bx = back;
				near.node = n.i0;
				far.bx = front;
				far.node = n.i1;
			}
			else {
				near.bx = front;
				near.node = n.i1;
				far.bx = back;
				far.node = n.i0;
			}
		}

		size_t next_dim = ((task.dim + 1) > PARAPROBE_ZAXIS) ? PARAPROBE_XAXIS : (task.dim + 1); //task.dim++;

		//assert( near.bx.outside_proximity(target) == 0.0 ); //##MK:::DEBUG
		near.d = 0.; //##MK:: should be not significant so set to ... ##### 0.;
		near.dim = next_dim;

		far.d = far.bx.outside_proximity(target);
		far.dim = next_dim;

		tasks[current_task] = far; //MK::near must be processed before far!
		++current_task;

		assert( static_cast<size_t>(current_task) < stack_size);
		tasks[current_task] = near;

	} while ( current_task != -1 );
	return best;
}


p3dm1 kd_tree::nearest_external(const p3dm1 target, const vector<p3dm1> & sortedpoints, apt_xyz epsball ) const
{
	const size_t stack_size = 32;
	const cuboid aabb(min, max);

	traverse_task tasks[stack_size] = {aabb, 0, 0, aabb.outside_proximity(target)}; //first box, first node, dim, distance to first
	int current_task = 0;
	apt_xyz currentbestd = epsball;

	p3dm1 best = p3dm1();
	do
	{
		const traverse_task &task = tasks[current_task];
		if (task.d > currentbestd) {
			--current_task;
			continue;
		}

		const node &n = nodes[task.node];
		if ( is_leaf(n) == true ) {
			//if n is a node n.i0 and n.i1 specify contiguous range of points on SORTED!! sortedpoints
			for(size_t i = n.i0; i <= n.i1; ++i) {
				//##MK::do SIMD
				const p3dm1 query = sortedpoints[i];
				const apt_xyz d = euclidean_sqrd(target, query);
				if ( d > currentbestd ) {
					continue;
				} else {
					//it cannot happen that in this query we find the same point again because the query operates on a different tree ie we assume that target IS not included in sortedpoints
					currentbestd = d;
					best = query;
				}
			}
			--current_task;
			continue;
		}

		traverse_task near;
		traverse_task far;
		//assert( n.dimdebug == task.dim );
		if ( task.dim == PARAPROBE_XAXIS ) {
			const apt_xyz splitposx = n.splitpos; //MK::using node-local piece of information rather than probing random memory position
			const cuboid left = cuboid( task.bx.min, p3d(splitposx, task.bx.max.y, task.bx.max.z) ); //#####MK::z coordinate of point...
			const cuboid right = cuboid( p3d(splitposx, task.bx.min.y, task.bx.min.z), task.bx.max );
			if ( target.x <= splitposx ) {
				near.bx = left;
				near.node = n.i0;
				far.bx = right;
				far.node = n.i1;
			}
			else {
				near.bx = right;
				near.node = n.i1;
				far.bx = left;
				far.node = n.i0;
			}
		}
		else if ( task.dim == PARAPROBE_YAXIS ) {
			const apt_xyz splitposy = n.splitpos;
			const cuboid lower = cuboid( task.bx.min, p3d(task.bx.max.x, splitposy, task.bx.max.z) );
			const cuboid upper = cuboid( p3d(task.bx.min.x, splitposy, task.bx.min.z), task.bx.max );
			if ( target.y <= splitposy ) {
				near.bx = lower;
				near.node = n.i0;
				far.bx = upper;
				far.node = n.i1;
			}
			else {
				near.bx = upper;
				near.node = n.i1;
				far.bx = lower;
				far.node = n.i0;
			}
		}
		else {
			const apt_xyz splitposz = n.splitpos;
			const cuboid back = cuboid( task.bx.min, p3d(task.bx.max.x, task.bx.max.y, splitposz) );
			const cuboid front = cuboid( p3d(task.bx.min.x, task.bx.min.y, splitposz), task.bx.max );
			if ( target.z <= splitposz ) {
				near.bx = back;
				near.node = n.i0;
				far.bx = front;
				far.node = n.i1;
			}
			else {
				near.bx = front;
				near.node = n.i1;
				far.bx = back;
				far.node = n.i0;
			}
		}

		size_t next_dim = ((task.dim + 1) > PARAPROBE_ZAXIS) ? PARAPROBE_XAXIS : (task.dim + 1); //task.dim++;

		//assert( near.bx.outside_proximity(target) == 0.0 ); //##MK:::DEBUG
		near.d = 0.; //##MK:: should be not significant so set to ... ##### 0.;
		near.dim = next_dim;

		far.d = far.bx.outside_proximity(target);
		far.dim = next_dim;

		tasks[current_task] = far; //MK::near must be processed before far!
		++current_task;

		assert( static_cast<size_t>(current_task) < stack_size);
		tasks[current_task] = near;

	} while ( current_task != -1 );
	return best;
}



//void kd_tree::range_rball_noclear_nosort(const vector<p3dm1> & sortedpoints, const apt_xyz radius_sqrd, vector<nbor> & result, const p3dm1 targ, size_t idxx = -1) //optional argument -1 if ##MK::for further nicefication of code make execute _external and local calls with the same underlying function
void kd_tree::range_rball_noclear_nosort_external(const p3dm1 target, const vector<p3dm1> & sortedpoints, const apt_xyz radius_sqrd, vector<nbor> & result )
{
	//MK::perform a spatial range query with a KDTree of a neighboring thread
	//MK::for these threadlocal KDTrees though may query points outside the neighboring threads local point portion, in that case
	//sortedpoints does not contain the target
	const size_t stack_size = 32;
	const cuboid aabb(min, max);
	//target is passed because not contained in sortedpoints of neighboring thread, const p3dm1 target = sortedpoints[idx];

	traverse_task tasks[stack_size] = {aabb, 0, 0, aabb.outside_proximity(target)}; //first box, first node, dim, distance to first
	int current_task = 0;

	//MK::function does neither clear results array beforehand nor sorts it afterwards!
	do
	{
		const traverse_task &task = tasks[current_task];
		if (task.d > radius_sqrd) {
			--current_task;
			continue;
		}

		const node &n = nodes[task.node];
		if ( is_leaf(n) == true ) {
			//if n is a node n.i0 and n.i1 specify contiguous range of points on SORTED!! sortedpoints
			for(size_t i = n.i0; i <= n.i1; ++i) {
				//##MK::do SIMD
				const p3dm1 query = sortedpoints[i];
				const apt_xyz d = euclidean_sqrd(target, query);
				if ( d > radius_sqrd ) {
					continue;
				} else {
					//if ( target != query ) {//MK::cannot happen that target is query because target is not in sortedpoints, therefore the function is external
						result.push_back( nbor(sqrt(d), query.m) );
					//}
				}
			}
			--current_task;
			continue;
		}

		traverse_task near;
		traverse_task far;
		//assert( n.dimdebug == task.dim );
		if ( task.dim == PARAPROBE_XAXIS ) {
			const apt_xyz splitposx = n.splitpos; //MK::using node-local piece of information rather than probing random memory position
			const cuboid left = cuboid( task.bx.min, p3d(splitposx, task.bx.max.y, task.bx.max.z) ); //#####MK::z cordinate of point...
			const cuboid right = cuboid( p3d(splitposx, task.bx.min.y, task.bx.min.z), task.bx.max );
			if ( target.x <= splitposx ) {
				near.bx = left;
				near.node = n.i0;
				far.bx = right;
				far.node = n.i1;
			}
			else {
				near.bx = right;
				near.node = n.i1;
				far.bx = left;
				far.node = n.i0;
			}
		}
		else if ( task.dim == PARAPROBE_YAXIS ) {
			const apt_xyz splitposy = n.splitpos;
			const cuboid lower = cuboid( task.bx.min, p3d(task.bx.max.x, splitposy, task.bx.max.z) );
			const cuboid upper = cuboid( p3d(task.bx.min.x, splitposy, task.bx.min.z), task.bx.max );
			if ( target.y <= splitposy ) {
				near.bx = lower;
				near.node = n.i0;
				far.bx = upper;
				far.node = n.i1;
			}
			else {
				near.bx = upper;
				near.node = n.i1;
				far.bx = lower;
				far.node = n.i0;
			}
		}
		else {
			const apt_xyz splitposz = n.splitpos;
			const cuboid back = cuboid( task.bx.min, p3d(task.bx.max.x, task.bx.max.y, splitposz) );
			const cuboid front = cuboid( p3d(task.bx.min.x, task.bx.min.y, splitposz), task.bx.max );
			if ( target.z <= splitposz ) {
				near.bx = back;
				near.node = n.i0;
				far.bx = front;
				far.node = n.i1;
			}
			else {
				near.bx = front;
				near.node = n.i1;
				far.bx = back;
				far.node = n.i0;
			}
		}

		size_t next_dim = ((task.dim + 1) > PARAPROBE_ZAXIS) ? PARAPROBE_XAXIS : (task.dim + 1); //task.dim++;

		//assert( near.bx.outside_proximity(target) == 0.0 ); //##MK:::DEBUG
		near.d = 0.; //##MK:: should be not significant so set to ... ##### 0.;
		near.dim = next_dim;

		far.d = far.bx.outside_proximity(target);
		far.dim = next_dim;

		tasks[current_task] = far; //MK::near must be processed before far!
		++current_task;

		assert( static_cast<size_t>(current_task) < stack_size);
		tasks[current_task] = near;

	} while ( current_task != -1 );
}




void kd_tree::range_rball_noclear_nosort(const size_t idx, const vector<p3dm1> & sortedpoints, const apt_xyz radius_sqrd, vector<nbor> & result )
{
	//MK::performs a spatial range query with a KDTree for which the target point is in the sortedpoints --> applies to domain global KDTrees
	//MK::threadlocal KDTrees though may query points outside the neighboring threads local point portion, in that case
	//sortedpoints does not contain the target
	const size_t stack_size = 32;
	const cuboid aabb(min, max);
	const p3dm1 target = sortedpoints[idx];

	traverse_task tasks[stack_size] = {aabb, 0, 0, aabb.outside_proximity(target)}; //first box, first node, dim, distance to first
	int current_task = 0;

	//MK::function does neither clear results array beforehand nor sorts it afterwards!
	do
	{
		const traverse_task &task = tasks[current_task];
		if (task.d > radius_sqrd) {
			--current_task;
			continue;
		}

		const node &n = nodes[task.node];
		if ( is_leaf(n) == true ) {
			//if n is a node n.i0 and n.i1 specify contiguous range of points on SORTED!! sortedpoints
			for(size_t i = n.i0; i <= n.i1; ++i) {
				//##MK::do SIMD
				const p3dm1 query = sortedpoints[i];
				const apt_xyz d = euclidean_sqrd(target, query);
				if ( d > radius_sqrd ) {
					continue;
				} else {
					if ( idx != i ) { //do not report myself!
						result.push_back( nbor(sqrt(d), query.m) );
					}
				}
			}
			--current_task;
			continue;
		}

		traverse_task near;
		traverse_task far;
		//assert( n.dimdebug == task.dim );
		if ( task.dim == PARAPROBE_XAXIS ) {
			const apt_xyz splitposx = n.splitpos; //MK::using node-local piece of information rather than probing random memory position
			const cuboid left = cuboid( task.bx.min, p3d(splitposx, task.bx.max.y, task.bx.max.z) ); //#####MK::z cordinate of point...
			const cuboid right = cuboid( p3d(splitposx, task.bx.min.y, task.bx.min.z), task.bx.max );
			if ( target.x <= splitposx ) {
				near.bx = left;
				near.node = n.i0;
				far.bx = right;
				far.node = n.i1;
			}
			else {
				near.bx = right;
				near.node = n.i1;
				far.bx = left;
				far.node = n.i0;
			}
		}
		else if ( task.dim == PARAPROBE_YAXIS ) {
			const apt_xyz splitposy = n.splitpos;
			const cuboid lower = cuboid( task.bx.min, p3d(task.bx.max.x, splitposy, task.bx.max.z) );
			const cuboid upper = cuboid( p3d(task.bx.min.x, splitposy, task.bx.min.z), task.bx.max );
			if ( target.y <= splitposy ) {
				near.bx = lower;
				near.node = n.i0;
				far.bx = upper;
				far.node = n.i1;
			}
			else {
				near.bx = upper;
				near.node = n.i1;
				far.bx = lower;
				far.node = n.i0;
			}
		}
		else {
			const apt_xyz splitposz = n.splitpos;
			const cuboid back = cuboid( task.bx.min, p3d(task.bx.max.x, task.bx.max.y, splitposz) );
			const cuboid front = cuboid( p3d(task.bx.min.x, task.bx.min.y, splitposz), task.bx.max );
			if ( target.z <= splitposz ) {
				near.bx = back;
				near.node = n.i0;
				far.bx = front;
				far.node = n.i1;
			}
			else {
				near.bx = front;
				near.node = n.i1;
				far.bx = back;
				far.node = n.i0;
			}
		}

		size_t next_dim = ((task.dim + 1) > PARAPROBE_ZAXIS) ? PARAPROBE_XAXIS : (task.dim + 1); //task.dim++;

		//assert( near.bx.outside_proximity(target) == 0.0 ); //##MK:::DEBUG
		near.d = 0.; //##MK:: should be not significant so set to ... ##### 0.;
		near.dim = next_dim;

		far.d = far.bx.outside_proximity(target);
		far.dim = next_dim;

		tasks[current_task] = far; //MK::near must be processed before far!
		++current_task;

		assert( static_cast<size_t>(current_task) < stack_size);
		tasks[current_task] = near;

	} while ( current_task != -1 );
}






void kd_tree::range_rball_noclear_nosort_external_p3dm1(const p3dm1 target, const vector<p3dm1> & sortedpoints, const apt_xyz radius_sqrd, vector<p3dm1> & result )
{
	//MK::perform a spatial range query with a KDTree of a neighboring thread
	//MK::for these threadlocal KDTrees though may query points outside the neighboring threads local point portion, in that case
	//sortedpoints does not contain the target
	const size_t stack_size = 32;
	const cuboid aabb(min, max);
	//target is passed because not contained in sortedpoints of neighboring thread, const p3dm1 target = sortedpoints[idx];

	traverse_task tasks[stack_size] = {aabb, 0, 0, aabb.outside_proximity(target)}; //first box, first node, dim, distance to first
	int current_task = 0;

	//MK::function does neither clear results array beforehand nor sorts it afterwards!
	do
	{
		const traverse_task &task = tasks[current_task];
		if (task.d > radius_sqrd) {
			--current_task;
			continue;
		}

		const node &n = nodes[task.node];
		if ( is_leaf(n) == true ) {
			//if n is a node n.i0 and n.i1 specify contiguous range of points on SORTED!! sortedpoints
			for(size_t i = n.i0; i <= n.i1; ++i) {
				//##MK::do SIMD
				const p3dm1 query = sortedpoints[i];
				const apt_xyz d = euclidean_sqrd(target, query);
				if ( d > radius_sqrd ) {
					continue;
				} else {
					//if ( target != query ) {//MK::cannot happen that target is query because target is not in sortedpoints, therefore the function is external
						result.push_back( query );
					//}
				}
			}
			--current_task;
			continue;
		}

		traverse_task near;
		traverse_task far;
		//assert( n.dimdebug == task.dim );
		if ( task.dim == PARAPROBE_XAXIS ) {
			const apt_xyz splitposx = n.splitpos; //MK::using node-local piece of information rather than probing random memory position
			const cuboid left = cuboid( task.bx.min, p3d(splitposx, task.bx.max.y, task.bx.max.z) ); //#####MK::z cordinate of point...
			const cuboid right = cuboid( p3d(splitposx, task.bx.min.y, task.bx.min.z), task.bx.max );
			if ( target.x <= splitposx ) {
				near.bx = left;
				near.node = n.i0;
				far.bx = right;
				far.node = n.i1;
			}
			else {
				near.bx = right;
				near.node = n.i1;
				far.bx = left;
				far.node = n.i0;
			}
		}
		else if ( task.dim == PARAPROBE_YAXIS ) {
			const apt_xyz splitposy = n.splitpos;
			const cuboid lower = cuboid( task.bx.min, p3d(task.bx.max.x, splitposy, task.bx.max.z) );
			const cuboid upper = cuboid( p3d(task.bx.min.x, splitposy, task.bx.min.z), task.bx.max );
			if ( target.y <= splitposy ) {
				near.bx = lower;
				near.node = n.i0;
				far.bx = upper;
				far.node = n.i1;
			}
			else {
				near.bx = upper;
				near.node = n.i1;
				far.bx = lower;
				far.node = n.i0;
			}
		}
		else {
			const apt_xyz splitposz = n.splitpos;
			const cuboid back = cuboid( task.bx.min, p3d(task.bx.max.x, task.bx.max.y, splitposz) );
			const cuboid front = cuboid( p3d(task.bx.min.x, task.bx.min.y, splitposz), task.bx.max );
			if ( target.z <= splitposz ) {
				near.bx = back;
				near.node = n.i0;
				far.bx = front;
				far.node = n.i1;
			}
			else {
				near.bx = front;
				near.node = n.i1;
				far.bx = back;
				far.node = n.i0;
			}
		}

		size_t next_dim = ((task.dim + 1) > PARAPROBE_ZAXIS) ? PARAPROBE_XAXIS : (task.dim + 1); //task.dim++;

		//assert( near.bx.outside_proximity(target) == 0.0 ); //##MK:::DEBUG
		near.d = 0.; //##MK:: should be not significant so set to ... ##### 0.;
		near.dim = next_dim;

		far.d = far.bx.outside_proximity(target);
		far.dim = next_dim;

		tasks[current_task] = far; //MK::near must be processed before far!
		++current_task;

		assert( static_cast<size_t>(current_task) < stack_size);
		tasks[current_task] = near;

	} while ( current_task != -1 );
}


void kd_tree::range_rball_noclear_nosort_p3dm1(const size_t idx, const vector<p3dm1> & sortedpoints, const apt_xyz radius_sqrd, vector<p3dm1> & result )
{
	//MK::performs a spatial range query with a KDTree for which the target point is in the sortedpoints --> applies to domain global KDTrees
	//MK::threadlocal KDTrees though may query points outside the neighboring threads local point portion, in that case
	//sortedpoints does not contain the target
	const size_t stack_size = 32;
	const cuboid aabb(min, max);
	const p3dm1 target = sortedpoints[idx];

	traverse_task tasks[stack_size] = {aabb, 0, 0, aabb.outside_proximity(target)}; //first box, first node, dim, distance to first
	int current_task = 0;

	//MK::function does neither clear results array beforehand nor sorts it afterwards!
	do
	{
		const traverse_task &task = tasks[current_task];
		if (task.d > radius_sqrd) {
			--current_task;
			continue;
		}

		const node &n = nodes[task.node];
		if ( is_leaf(n) == true ) {
			//if n is a node n.i0 and n.i1 specify contiguous range of points on SORTED!! sortedpoints
			for(size_t i = n.i0; i <= n.i1; ++i) {
				//##MK::do SIMD
				const p3dm1 query = sortedpoints[i];
				const apt_xyz d = euclidean_sqrd(target, query);
				if ( d > radius_sqrd ) {
					continue;
				} else {
					if ( idx != i ) { //do not report myself!
						result.push_back( query );
					}
				}
			}
			--current_task;
			continue;
		}

		traverse_task near;
		traverse_task far;
		//assert( n.dimdebug == task.dim );
		if ( task.dim == PARAPROBE_XAXIS ) {
			const apt_xyz splitposx = n.splitpos; //MK::using node-local piece of information rather than probing random memory position
			const cuboid left = cuboid( task.bx.min, p3d(splitposx, task.bx.max.y, task.bx.max.z) ); //#####MK::z cordinate of point...
			const cuboid right = cuboid( p3d(splitposx, task.bx.min.y, task.bx.min.z), task.bx.max );
			if ( target.x <= splitposx ) {
				near.bx = left;
				near.node = n.i0;
				far.bx = right;
				far.node = n.i1;
			}
			else {
				near.bx = right;
				near.node = n.i1;
				far.bx = left;
				far.node = n.i0;
			}
		}
		else if ( task.dim == PARAPROBE_YAXIS ) {
			const apt_xyz splitposy = n.splitpos;
			const cuboid lower = cuboid( task.bx.min, p3d(task.bx.max.x, splitposy, task.bx.max.z) );
			const cuboid upper = cuboid( p3d(task.bx.min.x, splitposy, task.bx.min.z), task.bx.max );
			if ( target.y <= splitposy ) {
				near.bx = lower;
				near.node = n.i0;
				far.bx = upper;
				far.node = n.i1;
			}
			else {
				near.bx = upper;
				near.node = n.i1;
				far.bx = lower;
				far.node = n.i0;
			}
		}
		else {
			const apt_xyz splitposz = n.splitpos;
			const cuboid back = cuboid( task.bx.min, p3d(task.bx.max.x, task.bx.max.y, splitposz) );
			const cuboid front = cuboid( p3d(task.bx.min.x, task.bx.min.y, splitposz), task.bx.max );
			if ( target.z <= splitposz ) {
				near.bx = back;
				near.node = n.i0;
				far.bx = front;
				far.node = n.i1;
			}
			else {
				near.bx = front;
				near.node = n.i1;
				far.bx = back;
				far.node = n.i0;
			}
		}

		size_t next_dim = ((task.dim + 1) > PARAPROBE_ZAXIS) ? PARAPROBE_XAXIS : (task.dim + 1); //task.dim++;

		//assert( near.bx.outside_proximity(target) == 0.0 ); //##MK:::DEBUG
		near.d = 0.; //##MK:: should be not significant so set to ... ##### 0.;
		near.dim = next_dim;

		far.d = far.bx.outside_proximity(target);
		far.dim = next_dim;

		tasks[current_task] = far; //MK::near must be processed before far!
		++current_task;

		assert( static_cast<size_t>(current_task) < stack_size);
		tasks[current_task] = near;

	} while ( current_task != -1 );
}


void kd_tree::range_rball_noclear_nosort_p3d( const p3dm1 target, vector<p3dm1> const & sortedpoints, const apt_xyz radius_sqrd, vector<p3dm1> & result )
{
	//MK::performs a spatial range query with a KDTree for which the target point is not in the sortedpoints
	const size_t stack_size = 32;
	const cuboid aabb(min, max);

	traverse_task tasks[stack_size] = {aabb, 0, 0, aabb.outside_proximity(target)}; //first box, first node, dim, distance to first
	int current_task = 0;

	//MK::function does neither clear results array beforehand nor sorts it afterwards!
	do
	{
		const traverse_task &task = tasks[current_task];
		if (task.d > radius_sqrd) {
			--current_task;
			continue;
		}

		const node &n = nodes[task.node];
		if ( is_leaf(n) == true ) {
			//if n is a node n.i0 and n.i1 specify contiguous range of points on SORTED!! sortedpoints
			for(size_t i = n.i0; i <= n.i1; ++i) {
				//##MK::do SIMD
				const p3dm1 query = sortedpoints[i];
				const apt_xyz d = euclidean_sqrd(target, query);
				if ( d > radius_sqrd ) {
					continue;
				} else {
					//cannot be myself!
					result.push_back( query );
				}
			}
			--current_task;
			continue;
		}

		traverse_task near;
		traverse_task far;
		//assert( n.dimdebug == task.dim );
		if ( task.dim == PARAPROBE_XAXIS ) {
			const apt_xyz splitposx = n.splitpos; //MK::using node-local piece of information rather than probing random memory position
			const cuboid left = cuboid( task.bx.min, p3d(splitposx, task.bx.max.y, task.bx.max.z) ); //#####MK::z cordinate of point...
			const cuboid right = cuboid( p3d(splitposx, task.bx.min.y, task.bx.min.z), task.bx.max );
			if ( target.x <= splitposx ) {
				near.bx = left;
				near.node = n.i0;
				far.bx = right;
				far.node = n.i1;
			}
			else {
				near.bx = right;
				near.node = n.i1;
				far.bx = left;
				far.node = n.i0;
			}
		}
		else if ( task.dim == PARAPROBE_YAXIS ) {
			const apt_xyz splitposy = n.splitpos;
			const cuboid lower = cuboid( task.bx.min, p3d(task.bx.max.x, splitposy, task.bx.max.z) );
			const cuboid upper = cuboid( p3d(task.bx.min.x, splitposy, task.bx.min.z), task.bx.max );
			if ( target.y <= splitposy ) {
				near.bx = lower;
				near.node = n.i0;
				far.bx = upper;
				far.node = n.i1;
			}
			else {
				near.bx = upper;
				near.node = n.i1;
				far.bx = lower;
				far.node = n.i0;
			}
		}
		else {
			const apt_xyz splitposz = n.splitpos;
			const cuboid back = cuboid( task.bx.min, p3d(task.bx.max.x, task.bx.max.y, splitposz) );
			const cuboid front = cuboid( p3d(task.bx.min.x, task.bx.min.y, splitposz), task.bx.max );
			if ( target.z <= splitposz ) {
				near.bx = back;
				near.node = n.i0;
				far.bx = front;
				far.node = n.i1;
			}
			else {
				near.bx = front;
				near.node = n.i1;
				far.bx = back;
				far.node = n.i0;
			}
		}

		size_t next_dim = ((task.dim + 1) > PARAPROBE_ZAXIS) ? PARAPROBE_XAXIS : (task.dim + 1); //task.dim++;

		//assert( near.bx.outside_proximity(target) == 0.0 ); //##MK:::DEBUG
		near.d = 0.; //##MK:: should be not significant so set to ... ##### 0.;
		near.dim = next_dim;

		far.d = far.bx.outside_proximity(target);
		far.dim = next_dim;

		tasks[current_task] = far; //MK::near must be processed before far!
		++current_task;

		assert( static_cast<size_t>(current_task) < stack_size);
		tasks[current_task] = near;

	} while ( current_task != -1 );
}







void kd_tree::range_rball_noclear_nosort_indices(const size_t idx, const vector<p3dm1> & sortedpoints, const apt_xyz radius_sqrd, vector<size_t> & result )
{
	//MK::performs a spatial range query with a KDTree for which the target point is in the sortedpoints --> applies to domain global KDTrees
	//MK::threadlocal KDTrees though may query points outside the neighboring threads local point portion, in that case
	//sortedpoints does not contain the target
	const size_t stack_size = 32;
	const cuboid aabb(min, max);
	const p3dm1 target = sortedpoints[idx];

	traverse_task tasks[stack_size] = {aabb, 0, 0, aabb.outside_proximity(target)}; //first box, first node, dim, distance to first
	int current_task = 0;

	//MK::function does neither clear results array beforehand nor sorts it afterwards!
	do
	{
		const traverse_task &task = tasks[current_task];
		if (task.d > radius_sqrd) {
			--current_task;
			continue;
		}

		const node &n = nodes[task.node];
		if ( is_leaf(n) == true ) {
			//if n is a node n.i0 and n.i1 specify contiguous range of points on SORTED!! sortedpoints
			for(size_t i = n.i0; i <= n.i1; ++i) {
				//##MK::do SIMD
				const p3dm1 query = sortedpoints[i];
				const apt_xyz d = euclidean_sqrd(target, query);
				if ( d > radius_sqrd ) {
					continue;
				} else {
					if ( idx != i ) { //do not report myself!
						result.push_back( i );
					}
				}
			}
			--current_task;
			continue;
		}

		traverse_task near;
		traverse_task far;
		//assert( n.dimdebug == task.dim );
		if ( task.dim == PARAPROBE_XAXIS ) {
			const apt_xyz splitposx = n.splitpos; //MK::using node-local piece of information rather than probing random memory position
			const cuboid left = cuboid( task.bx.min, p3d(splitposx, task.bx.max.y, task.bx.max.z) ); //#####MK::z cordinate of point...
			const cuboid right = cuboid( p3d(splitposx, task.bx.min.y, task.bx.min.z), task.bx.max );
			if ( target.x <= splitposx ) {
				near.bx = left;
				near.node = n.i0;
				far.bx = right;
				far.node = n.i1;
			}
			else {
				near.bx = right;
				near.node = n.i1;
				far.bx = left;
				far.node = n.i0;
			}
		}
		else if ( task.dim == PARAPROBE_YAXIS ) {
			const apt_xyz splitposy = n.splitpos;
			const cuboid lower = cuboid( task.bx.min, p3d(task.bx.max.x, splitposy, task.bx.max.z) );
			const cuboid upper = cuboid( p3d(task.bx.min.x, splitposy, task.bx.min.z), task.bx.max );
			if ( target.y <= splitposy ) {
				near.bx = lower;
				near.node = n.i0;
				far.bx = upper;
				far.node = n.i1;
			}
			else {
				near.bx = upper;
				near.node = n.i1;
				far.bx = lower;
				far.node = n.i0;
			}
		}
		else {
			const apt_xyz splitposz = n.splitpos;
			const cuboid back = cuboid( task.bx.min, p3d(task.bx.max.x, task.bx.max.y, splitposz) );
			const cuboid front = cuboid( p3d(task.bx.min.x, task.bx.min.y, splitposz), task.bx.max );
			if ( target.z <= splitposz ) {
				near.bx = back;
				near.node = n.i0;
				far.bx = front;
				far.node = n.i1;
			}
			else {
				near.bx = front;
				near.node = n.i1;
				far.bx = back;
				far.node = n.i0;
			}
		}

		size_t next_dim = ((task.dim + 1) > PARAPROBE_ZAXIS) ? PARAPROBE_XAXIS : (task.dim + 1); //task.dim++;

		//assert( near.bx.outside_proximity(target) == 0.0 ); //##MK:::DEBUG
		near.d = 0.; //##MK:: should be not significant so set to ... ##### 0.;
		near.dim = next_dim;

		far.d = far.bx.outside_proximity(target);
		far.dim = next_dim;

		tasks[current_task] = far; //MK::near must be processed before far!
		++current_task;

		assert( static_cast<size_t>(current_task) < stack_size);
		tasks[current_task] = near;

	} while ( current_task != -1 );
}


inline p3d kd_tree::get_min()
{
	return this->min;
}


inline p3d kd_tree::get_max()
{
	return this->max;
}


bool kd_tree::verify(const vector<p3dm1> & sortedpoints)
{
	vector<unsigned char> visited;
	visited.reserve(sortedpoints.size());
	for(size_t i = 0; i < sortedpoints.size(); ++i) { visited.push_back(0x00); }
	for(size_t i = 0; i < nodes.size(); ++i) { //check that point references for each node do not overlap!
		node &n = nodes[i];
		if(is_leaf(n) == true) {
			for(size_t k = n.i0; k <= n.i1; ++k) { //n.i0 and n.i1 are inclusive!
				if(visited[k] == 0x00)
					visited[k] = 0x01;
				else
					return false;
			}
		}
	}
	return true;
}


void kd_tree::get_allboundingboxes( vector<scan_task> & out )
{
	//returns all bounding boxes of the tree
	const size_t stack_size = 32;
	cuboid aabb(min, max);
	scan_task tasks[stack_size] = {aabb, 0, 0}; //first box, first node, dim, distance to first
	int current_task = 0;
	out.push_back( tasks[current_task] );

	do
	{
		const scan_task &task = tasks[current_task];
		const node &n = nodes[task.node];

		if ( is_leaf(n) == true ) {
			//if arrived at leaf we just have to jump back, the corresponding box pair was written out already at one level higher before
			--current_task;
			continue;
		}

		scan_task aa; //##MK::irrelevant whether aa or bb are qualitatively near or far we want to visit all nodes...
		scan_task bb;
		if ( task.dim == PARAPROBE_XAXIS ) {
			const apt_xyz splitposx = n.splitpos;
			const cuboid a = cuboid( task.bx.min, p3d(splitposx, task.bx.max.y, task.bx.max.z) );
			const cuboid b = cuboid( p3d(splitposx, task.bx.min.y, task.bx.min.z), task.bx.max );
			aa.bx = a;
			aa.node = n.i0;
			bb.bx = b;
			bb.node = n.i1;
		}
		else if ( task.dim == PARAPROBE_YAXIS ) {
			const apt_xyz splitposy = n.splitpos;
			const cuboid a = cuboid( task.bx.min, p3d(task.bx.max.x, splitposy, task.bx.max.z) );
			const cuboid b = cuboid( p3d(task.bx.min.x, splitposy, task.bx.min.z), task.bx.max );
			aa.bx = a;
			aa.node = n.i0;
			bb.bx = b;
			bb.node = n.i1;
		}
		else {
			const apt_xyz splitposz = n.splitpos;
			const cuboid a = cuboid( task.bx.min, p3d(task.bx.max.x, task.bx.max.y, splitposz) );
			const cuboid b = cuboid( p3d(task.bx.min.x, task.bx.min.y, splitposz), task.bx.max );
			aa.bx = a;
			aa.node = n.i0;
			bb.bx = b;
			bb.node = n.i1;
		}

		size_t next_dim = ((task.dim + 1) > PARAPROBE_ZAXIS) ? PARAPROBE_XAXIS : (task.dim + 1); //circular looping
		aa.dim = next_dim;
		bb.dim = next_dim;

		tasks[current_task] = aa; //MK::near must be processed before far!
		out.push_back( aa );
		++current_task;

		assert( static_cast<size_t>(current_task) < stack_size);
		out.push_back( bb );
		tasks[current_task] = bb;

	} while ( current_task != -1 );
}


void kd_tree::display_nodes()
{
	for(size_t i = 0; i < nodes.size(); ++i) {
		node &n = nodes[i];
		if(is_leaf(n)==true) {
			cout << n.i0 << "\t\t" << n.i1 << "\t\t" << n.split << "\t\t" << (nodes[i].i1 - nodes[i].i0 + 1) << "\n";
		}
	}
}


size_t kd_tree::get_treememory_consumption()
{
	size_t bytes = 2*sizeof(p3d) + nodes.size()*sizeof(node);
	return bytes; //static_cast<double>(bytes) / (1024 * 1024); //MB
}
