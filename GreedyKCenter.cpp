//http://www.mathworks.com/help/matlab/matlab_external/debugging-c-c-language-mex-files.html

//Programmer: Chris Tralie
//Purpose: To implement the naive O(lambda^(O(1)) n log Phi) algorithm to perform the
//greedy 2-approximation to the k center ordering
#include <mex.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <iostream>
#include <algorithm>
#include <queue>
#include <list>
#include <vector>
#include <assert.h>

using namespace std;

#define VERBOSE 1
#define VERBOSE_PRINTCENTERS 1
#define VERBOSE_HEAP 0

//TODO: Accelerate with GPU
//TODO: Add a distance oracle in place of requiring things to be Euclidean


//NOTE: As suggested in the paper, points from friends lists are removed 
//in a lazy fashion whenever the lists are scanned.


//Structure used in the max heap to find the furthest point from the 
//k centers, which will be used as the next center in the greedy scheme
class FPElem {
public:
	int k; //Associated with this cluster center
	int index; //Index of the furthest point in the points list
	double dist; //Distance of the furthest point
	
	void print() {
		mexPrintf("k = %i, index = %i, dist = %g\n", k, index, dist);
	}
};

//Comparator for this structure
class FPComp {
public:
	bool operator()(const FPElem& e1, const FPElem& e2) {
		if (e1.dist < e2.dist) {
			return true;
		}
		return false;
	}
};


//Class for storing each k-center along with its friends list, cluster list
//and some other auxiliary information
class KCenter {
public:
	int k;
	int pointIndex;
	double R;
	int furthestP;//Index of furthest point in cluster
	double furthestD;//Distance of furthest point in cluster
	list<int> points;//Points that belong to this cluster (indexed by occurrence in X)
	list<int> friends;//Friend cluster centers (indexed by KCenter.k)
	
	KCenter() {
		k = 0;
		pointIndex = 0;
		R = 0;
		furthestP = 0;
		furthestD = 0;
	}
};

//Holds all of the data structures that are used as the problem evolves
class ProblemInstance {
public:
	double* alphasSqr;//Squared distance of points to cluster centers
	KCenter* centers;//The k centers
	priority_queue<FPElem, vector<FPElem>, FPComp> FPHeap;//The structure for containing the max heap for distances
	double* X;//The actual N-point point cloud in R^D
	size_t N, D;
	
	int phase;//Which phase the algorithm is in
	double phaseR;//The radius of the phase the algorithm is in
	//The cluster centers which serve the points at the end of phases
	int* lastServingCenter;
	int* lastLastServingCenter;
	
	ProblemInstance(double* paramX, size_t paramN, size_t paramD) {
		X = paramX;
		N = paramN;
		D = paramD;
		alphasSqr = new double[N];
		centers =  new KCenter[N];
		phase = 0;
		phaseR = 0;
		lastServingCenter = new int[N];
		lastLastServingCenter = new int[N];
	}
	
	~ProblemInstance() {
		delete[] alphasSqr;
		delete[] centers;
		delete[] lastServingCenter;
		delete[] lastLastServingCenter;
	}
	
	void copyCentersRads(double* c, double* r) {
		for (size_t i = 0; i < N; i++) {
			c[i] = centers[i].pointIndex + 1;//Matlab is 1-indexed
			r[i] = centers[i].R;
		}
	}
	
	void printCenter(size_t k) {
		list<int>::iterator it;
		mexPrintf("===== Cluster %i =====\n", centers[k].pointIndex);
		mexPrintf("Radius %g, furthestP %i, furthestD %g\n", centers[k].R, centers[k].furthestP, centers[k].furthestD);
		mexPrintf("Members: ");
		it = centers[k].points.begin();
		while (it != centers[k].points.end()) {
			mexPrintf("%i, ", *it);
			it++;
		}
		mexPrintf("\nFriends: ");
		it = centers[k].friends.begin();
		while (it != centers[k].friends.end()) {
			mexPrintf("%i, ", *it);
			it++;
		}
		mexPrintf("\n\n");		
	}
	
	//For debugging
	void printAllCenters(size_t k) {
		for (size_t i = 0; i <= k; i++) {
			printCenter(i);
		}
	}
};

//Return the squared distance between X_i and X_j
double getSqrDist(ProblemInstance& inst, int i, int j) {
	double ret = 0.0;
	double diff = 0.0;
	for (size_t k = 0; k < inst.D; k++) {
		diff = inst.X[k*inst.N + i] - inst.X[k*inst.N + j];
		ret += diff*diff;
	}
	return ret;
}

//Helper function that updates friends lists and served points of a cluster
//center clusterK and also moves points to the new cluster center newclusterK
//This function is called for every friend in the newclusterK's previous 
//center friends list
//Inputs: *inst: Problem instance
//*newclusterK: The index of the new cluster center into insts's cluster list
//*clusterK: The index of another cluster center that is being checked
void scanCluster(ProblemInstance& inst, int newclusterK, int clusterK) {
	//Step 1: Check all of the points in the cluster center to see if
	//they are closer to the new center or not
	bool furthestInvalid = false;
	list<int>::iterator it = inst.centers[clusterK].points.begin();
	while (it != inst.centers[clusterK].points.end()) {
		int P = *it;
		//Compute the distance of the point in this cluster to the new k-center
		double distSqr = getSqrDist(inst, P, inst.centers[newclusterK].pointIndex);
		//If the distance is less than the distance to the current cluster center
		//then this point moves to the new k-center's cluster
		if (distSqr < inst.alphasSqr[P]) {
			//Remove this point from the old center
			it = inst.centers[clusterK].points.erase(it);
			inst.alphasSqr[P] = distSqr;
			//Add this point to the new center if it isn't the new cluster center
			if (P != inst.centers[newclusterK].pointIndex) {
				inst.centers[newclusterK].points.push_back(P);
			}
			//The furthest point in this cluster is no longer
			//part of the cluster
			if (P == inst.centers[clusterK].furthestP) {
				furthestInvalid = true;
			}
		}
		else {
			it++;//Only move the list iterator by one if the previous
			//point was not removed from the list
		}
	}
	
	//Step 2: If the furthest point in the kth cluster has been removed,
	//figure out what the the new furthest point in the kth cluster is, 
	//and add that point to the max heap
	if (furthestInvalid && inst.centers[clusterK].points.size() > 0) {
		FPElem newElem;
		newElem.k = clusterK;
		newElem.index = 0;
		newElem.dist = 0.0;
		it = inst.centers[clusterK].points.begin();
		while (it != inst.centers[clusterK].points.end()) {
			int P = *it;
			if (inst.alphasSqr[P] > newElem.dist) {
				newElem.dist = inst.alphasSqr[P];
				newElem.index = P;
			}
			it++;
		}
		newElem.dist = sqrt(newElem.dist);
		inst.FPHeap.push(newElem);
		inst.centers[clusterK].furthestP = newElem.index;
		inst.centers[clusterK].furthestD = newElem.dist;
		if (VERBOSE_HEAP) {
			mexPrintf("Pushing: ");
			newElem.print();
		}
	}
}

//This function's job is to get the next cluster center and
//to update all of the friends lists and clusters of centers
//that are nearby to the new cluster
void getNextClusterCenter(ProblemInstance& inst, int k) {
	if (VERBOSE)
		mexPrintf("\n\n\nCalculating cluster center %i\n", k);
	
	list<int>::iterator it;
	inst.centers[k].k = k;
	//Step 1: Extract the maximum value from the heap
	//and set the next cluster center to be this point
	FPElem e;
	do {
		//Loop through while outdated furthest points are on the heap
		e = inst.FPHeap.top();
		inst.FPHeap.pop();
		if (VERBOSE_HEAP) {
			mexPrintf("Popping: ");
			e.print();
			mexPrintf("inst.centers[%i].furthestP = %i\n", e.k, inst.centers[e.k].furthestP);
		}
	}
	while (e.index != inst.centers[e.k].furthestP);
	
	//Set the next cluster center to be this point
	inst.centers[k].pointIndex = e.index;
	//Set the radius of the previous cluster center to be
	//this distance
	inst.centers[k-1].R = e.dist;
	int CPk = e.k;
	if (VERBOSE)
		mexPrintf("Cluster center %i is point %i, radius %g\n", k, e.index + 1, e.dist);
	//Add this new cluster center to the previous cluster center's friends list
	inst.centers[k].friends.push_back(CPk);
	inst.centers[CPk].friends.push_back(k);

	//Step 2: Scan all of the points currently served by the same
	//cluster as CPk and by clusters in CPk's friends list
	scanCluster(inst, k, CPk);
	it = inst.centers[CPk].friends.begin();
	while (it != inst.centers[CPk].friends.end()) {
		scanCluster(inst, k, inst.centers[*it].k);
		//Lazy check to see if this friend should be removed from the friends list
		double friendDist = sqrt(getSqrDist(inst, inst.centers[CPk].pointIndex, inst.centers[*it].pointIndex));
		if (friendDist > 8*e.dist) {
			mexPrintf("Removing %i from %i's friends list\n", *it, CPk);
			it = inst.centers[CPk].friends.erase(it);
		}
		else {
			it++;
		}
	}

	//Step 3: Determine the furthest point in the new cluster and add
	//it to the heap
	if (inst.centers[k].points.size() > 0) {
		it = inst.centers[k].points.begin();
		while (it != inst.centers[k].points.end()) {
			double distSqr = getSqrDist(inst, *it, inst.centers[k].pointIndex);
			if (distSqr > inst.centers[k].furthestD) {
				inst.centers[k].furthestD = distSqr;
				inst.centers[k].furthestP = *it;
			}
			it++;
		}
		inst.centers[k].furthestD = sqrt(inst.centers[k].furthestD);
		FPElem newClusterFar;
		newClusterFar.k = k;
		newClusterFar.index = inst.centers[k].furthestP;
		newClusterFar.dist = inst.centers[k].furthestD;
		inst.FPHeap.push(newClusterFar);
		if (VERBOSE_HEAP) {
			mexPrintf("Pushing: ");
			newClusterFar.print();
		}
	}


	//Step 4: Construct the friends list of the new center
	int Pk = e.index;//New center point index
	double Rk = e.dist;//(This is an upper bound of what was said in the paper)
	int LSk = inst.lastLastServingCenter[Pk];//Last serving center
	it = inst.centers[LSk].friends.begin();
	//Look at the friends list of the new cluster center 2 phases ago
	//and add all of the points that are at most 4*Rk from the new center
	//Also add the new center to those lists
	while (it != inst.centers[LSk].friends.end()) {
		//Distance of a friend to the serving center two phases ago
		double distCenter = sqrt(getSqrDist(inst, inst.centers[*it].pointIndex, LSk));
		//Distance to new cluster center
		double distNew = sqrt(getSqrDist(inst, inst.centers[*it].pointIndex, Pk));
		if (distNew < 4*Rk && *it != k) {
			//Add the centers to each others' friends lists
			if (VERBOSE)
				mexPrintf("Adding friends %i and %i\n", *it, k);
			inst.centers[k].friends.push_back(*it);
			inst.centers[*it].friends.push_back(k);
		}
		//Lazy check if this friend should be removed from the friends list
		if (distCenter > 8*Rk) {
			if (VERBOSE)
				mexPrintf("Removing %i from %i's friends list\n", *it, CPk);
			it = inst.centers[CPk].friends.erase(it);
		}
		else {
			it++;
		}		
	}

	//Step 5: See if it's time to move onto a new phase
	//(this step takes (O(n)) time and is hit O(log(Spread)) times)
	if (e.dist <= inst.phaseR/2) {
		inst.phase++;
		inst.phaseR = e.dist;
		if (VERBOSE)
			mexPrintf("Moving onto phase %i, Radius %g\n", inst.phase, inst.phaseR);
		//Update lastServingCenter and lastLastServingCenter
		for (size_t i = 0; i < inst.N; i++) {
			inst.lastLastServingCenter[i] = inst.lastServingCenter[i];
		}
		size_t NTotalCluster = 0;
		for (int i = 0; i <= k; i++) {
			NTotalCluster += inst.centers[i].points.size() + 1;
			it = inst.centers[i].points.begin();
			while (it != inst.centers[i].points.end()) {
				inst.lastServingCenter[*it] = i;
				it++;
			}
		}
		//Print error if the sum of the number of points belonging to
		//each cluster center does not equal N
		if (NTotalCluster != inst.N) {
			mexPrintf("Error: NTotalCluster (%i) != N (%i)\n", NTotalCluster, inst.N);
		}
	}
}

//Inputs: X (NxD matrix of points, where D is the dimension)
//Outputs: centers (1xN ordering of centers)
//rads (1xN matrix of radii needed to cover points at each level)
void mexFunction(int nOutArray, mxArray *OutArray[], int nInArray, const mxArray *InArray[]) {  
	///////////////MEX INPUTS/////////////////
	if (nInArray < 1) {
		mexErrMsgTxt("Distance Matrix Input Required");
		return;
	}
	
	const mwSize *dims;
	size_t N, D;
	int ndim = mxGetNumberOfDimensions(InArray[0]);
	if (ndim != 2) {
		mexErrMsgTxt("Expecting 2D matrix as first input");
		return;
	}
	dims = mxGetDimensions(InArray[0]);
	N = dims[0];
	D = dims[1];
	
	///////////////ALGORITHM/////////////////
	
	//Initialize all variables
	ProblemInstance  inst((double*)mxGetPr(InArray[0]), N, D);
	
	//Set the first cluster center to be the first point
	inst.centers[0].pointIndex = 0;
	inst.centers[0].k = 0;
	//Initialize the points that belong to this cluster to be all points
	//and in the process figure out which point is the furthest
	for (size_t i = 1; i < N; i++) {
		double dSqr = getSqrDist(inst, 0, i);
		inst.alphasSqr[i] = dSqr;
		inst.centers[0].points.push_back(i);
		if (dSqr > inst.centers[0].furthestD) {
			inst.centers[0].furthestD = dSqr;
			inst.centers[0].furthestP = i;
		}
	}
	inst.centers[0].furthestD = sqrt(inst.centers[0].furthestD);
	FPElem e0;
	e0.k = 0; 
	e0.index = inst.centers[0].furthestP; 
	e0.dist = inst.centers[0].furthestD;
	inst.FPHeap.push(e0);
	
	//Initialize phase information for the first phase
	inst.phase = 0;
	inst.phaseR = e0.dist;
	for (size_t i = 0; i < N; i++) {
		inst.lastServingCenter[i] = 0;
		inst.lastLastServingCenter[i] = 0;
	}
	
	//Loop through and compute the new cluster centers
	//one at a time
	for (size_t k = 1; k < N; k++) {
		if (VERBOSE_PRINTCENTERS)
			mexPrintf("*********STEP %i************\n", k);
		getNextClusterCenter(inst, k);
		if (VERBOSE_PRINTCENTERS)
			inst.printAllCenters(k);
	}
	
	///////////////MEX OUTPUTS/////////////////
	double* centers = new double[N];//Index of the centers chosen in the order they are chosen
	double* rads = new double[N];//Radius of the balls needed around the (i <= k) center needed
	//to cover all points at level k
	inst.copyCentersRads(centers, rads);
	
	mwSize outdims[2];
	outdims[0] = N;
	outdims[1] = 1;
	OutArray[0] = mxCreateNumericArray(2, outdims, mxDOUBLE_CLASS, mxREAL);
	double* centersOut = (double*)mxGetPr(OutArray[0]);
	OutArray[1] = mxCreateNumericArray(2, outdims, mxDOUBLE_CLASS, mxREAL);
	double* radsOut = (double*)mxGetPr(OutArray[1]);
	memcpy(centersOut, centers, N*sizeof(double));
	memcpy(radsOut, rads, N*sizeof(double));
	
	///////////////CLEANUP/////////////////
	delete[] centers;
	delete[] rads;
}
