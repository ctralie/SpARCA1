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
#include <vector>
#include <assert.h>

using namespace std;

//TODO: *Experiment with using square distances throughout on running time
//TODO: Accelerate with GPU
//TODO: See if vectors are slowing me down


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
	vector<int> points;//Points that belong to this cluster (indexed by occurrence in X)
	vector<int> friends;//Friend cluster centers (indexed by KCenter.k)
	
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
	double* alphasSqr;//Distance of points to cluster centers
	KCenter* centers;//The k centers
	vector<FPElem> FPHeap;//The structure for containing the min heap for distances
	double* X;
	size_t N, D;
	int phase;
	
	ProblemInstance(double* paramX, size_t paramN, size_t paramD) {
		X = paramX;
		N = paramN;
		D = paramD;
		alphasSqr = new double[N];
		centers =  new KCenter[N];
		phase = 0;
	}
	
	~ProblemInstance() {
		delete[] alphasSqr;
		delete[] centers;
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

//Helper function that updates friends lists and served points of a cluster center
//and also moves points to the new cluster center
void scanCluster(ProblemInstance& inst, int newclusterK, int clusterK) {
	vector<int> pointsToRemove;
	bool furthestInvalid = false;
	for (size_t i = 0; i < centers[clusterK].points.size(); i++) {
		int P = centers[clusterK].points[i];
		//Compute the distance of the point in this cluster to the new k-center
		double distSqr = getSqrDistance(inst.X, inst.N, inst.D, P, inst.centers[newclusterK].pointIndex);
		//If the distance is less than the distance to the current cluster center
		//then this point moves to the new k-center's cluster
		if (distSqr < alphasSqr[P]) {
			pointsToRemove.push_back(i);
			centers[newClusterK].points.push_back(P);
			alphasSqr[P] = distSqr;
			if (P == centers[clusterK].furthestP) {
				furthestInvalid = true;
				//The furthest point in this cluster is no longer
				//part of the cluster
			}
		}
	}
	//Now remove all of the points from the cluster list that are no
	//longer part of the cluster and move them to the new one
	for (int i = (int)pointsToRemove.size() - 1; i >= 0; i--) {
		int P = pointsToRemove[i];
		centers[clusterK].
	} 
}

void getNextClusterCenter(ProblemInstance& ins, int k) {
	inst.centers[k].k = k;
	//Step 1: Extract the maximum value from the heap
	//and set the next cluster center to be this point
	FPElem e;
	do {
		pop_heap(inst.FPHeap.begin(), inst.FPHeap.end(), &FPComp);
		e = inst.FPComp.pop_back();
		//Loop through while outdated closest points are on the heap
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
	
	//Loop through and compute the new cluster centers
	//one at a time
	for (size_t k = 1; k < N; k++) {
		getNextClusterCenter(inst, k);
	}
	
/*	///////////////MEX OUTPUTS/////////////////
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
