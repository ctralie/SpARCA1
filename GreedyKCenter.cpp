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
#include <list>
#include <assert.h>

using namespace std;

//TODO: Accelerate with GPU
//TODO: Add a distance oracle in place of requiring things to be Euclidean


//NOTE: As suggested in the paper, points from friends lists are removed 
//in a lazy fashion whenever the lists are scanned.

//Return the squared distance between X_i and X_j
double getSqrDist(double* X, size_t N, size_t D, int i, int j) {
	double ret = 0.0;
	double diff = 0.0;
	for (size_t k = 0; k < D; k++) {
		diff = D[k*N + i] - D[k*N + j];
		ret += diff*diff;
	}
	return ret;
}


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
	vector<FPElem> FPHeap;//The structure for containing the max heap for distances
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
}

//Structure used in the max heap to find the furthest point from the 
//k centers, which will be used as the next center in the greedy scheme
struct FPElem {
	int k; //Associated with this cluster center
	int index; //Index of the furthest point in the points list
	double dist; //Distance of the furthest point
};

//Comparator for this structure
//Sort in descending order
bool FPComp(const FPElem& e1, const FPElem& e2) {
	if (e2.dist < e1.dist) {
		return true;
	}
	return false;
}

//Helper function that updates friends lists and served points of a cluster
//center clusterK and also moves points to the new cluster center newclusterK

//This function is called for every friend in the newclusterK's previous 
//center friends list
void scanCluster(ProblemInstance& inst, int newclusterK, int clusterK) {
	bool furthestInvalid = false;
	vector<int>::iterator it = inst.centers[clusterK].points.begin();
	
	//Step 1: Check all of the points in the cluster center to see if
	//they are closer to the new center or not
	while (it != inst.centers[clusterK].points.end()) {
		int P = *it;
		//Compute the distance of the point in this cluster to the new k-center
		double distSqr = getSqrDistance(inst.X, inst.N, inst.D, P, inst.centers[newclusterK].pointIndex);
		//If the distance is less than the distance to the current cluster center
		//then this point moves to the new k-center's cluster
		if (distSqr < inst.alphasSqr[P]) {
			//Remove this point from the old center
			it = inst.centers[clusterK].points.erase(it);
			//Add this point to the new center
			inst.centers[newClusterK].points.push_back(P);
			inst.alphasSqr[P] = distSqr;
			if (P == inst.centers[clusterK].furthestP) {
				//The furthest point in this cluster is no longer
				//part of the cluster
				furthestInvalid = true;
			}
		}
		else {
			it++;//Only move the list iterator by one if the previous
			//point was not removed from the list
		}
	}
	
	//Step 2: If the furthest point in the kith cluster has been removed,
	//figure out what the the new furthest point in the kth cluster is, 
	//and add that point to the max heap
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
	}
	newElem.dist = sqrt(newElem.dist);
	inst.FPHeap.push_back(newElem);
	push_heap(inst.FPHeap);
}

//This function's job is to get the next cluster center and
//to update all of the friends lists and clusters of centers
//that are nearby to the new cluster
void getNextClusterCenter(ProblemInstance& ins, int k) {
	inst.centers[k].k = k;
	//Step 1: Extract the maximum value from the heap
	//and set the next cluster center to be this point
	FPElem e;
	do {
		//Loop through while outdated furthest points are on the heap
		pop_heap(inst.FPHeap.begin(), inst.FPHeap.end(), &FPComp);
		e = inst.FPComp.pop_back();
	}
	while (e.index != inst.centers[e.k].furthestP);
	//Set the next cluster center to be this point
	inst.centers[k].pointIndex = e.index;
	//Set the radius of the previous cluster center to be
	//this distance
	inst.centers[k-1].R = e.dist;
	int CPk = e.k;
	
	//Step 2: Scan all of the points currently served by the same
	//cluster as CPk and by clusters in CPk's friends list
	scanCluster(inst, k, CPk);
	for (size_t i = 0; i < inst.centers[CPk].friends.size(); i++) {
		scanCluster(inst, inst.centers[centers[CPk].friends[i]].k)
	}
	
	//Step 3: Construct the friends list of the new center
	//TODO:
	
	//Step 4: See if it's time to move onto a new phase
	//TODO:
	//(this step takes (O(n)) time and is hit O(log(Spread)) times)
	if (e.dist <= inst.phaseR/2) {
		inst.phase++;
		inst.phaseR = e.dist;
		
		//Print error if the sum of the number of points belonging to
		//each cluster center does not equal N
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
	for (size_t i = 0; i < N; i++) {
		float dSqr = getSqrDist(inst.X, N, D, 0, i);
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
	inst.FPHeap.push_back(e0);
	
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
		getNextClusterCenter(inst, k);
	}
	
/*	///////////////MEX OUTPUTS/////////////////
	//TODO
	double* centers = new double[N];//Index of the centers chosen in the order they are chosen
	double* rads = new double[N];//Radius of the balls needed around the (i <= k) center needed
	//to cover all points at level k
	
	mwSize outdims[2];
	//Output 0D homology classes
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
	delete[] rads;*/
	
}
