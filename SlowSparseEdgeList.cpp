//Programmer: Chris Tralie
//Purpose: To extract a sparse edge list from a point cloud given a cover tree
//over that point cloud
#include <mex.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <iostream>
#include <algorithm>
#include <vector>
#include <assert.h>
#include <map>


using namespace std;

//Return the squared distance between X_i and X_j
double getSqrDist(double* X, int N, int D, int i, int j) {
	double ret = 0.0;
	double diff = 0.0;
	for (size_t k = 0; k < D; k++) {
		diff = X[k*N + i] - X[k*N + j];
		ret += diff*diff;
	}
	return ret;
}

//Return by value e1, e2, ed.  Set e1 and e2 to -1 if the edge shouldn't be added
//This function figures out if and when an edge between vp and vq should be included
void getEdge(int p, int q, double* X, int N, int D, double* ts, int* e1, int* e2, double* ed) {
	double dpq = sqrt(getSqrDist(X, N, D, p, q));
	double tp = ts[p];
	double tq = ts[q];
	double alpha;
	
	if (tq < tp) {
		//Make sure tp indexes the earlier of the two death times
		double temp = tq;
		tq = tp;
		tp = temp;
	}
	//Case a
	if (dpq < tp/3) {
		*e1 = p;
		*e2 = q;
		*ed = dpq;
		//mexPrintf("e1 = %i, e2 = %i, dpq = %g, tp/3 = %g\n", p, q, dpq, tp/3);
		return;
	}
	//Case b
	alpha = 2*dpq - tp/3;
	if (alpha >= tp/3 && alpha <= tp && alpha <= tq/3) {
		*e1 = p;
		*e2 = q;
		*ed = alpha;
		//mexPrintf("e1 = %i, e2 = %i, alpha = %g, tp/3 = %g, tq/3 = %g\n", p, q, alpha, tp/3, tq/3);
		return;
	}
	
	//Otherwise, it's never added
	*e1 = -1;
	*e2 = -1;
}

//Inputs: *X: NxD point cloud matrix
//*radii: L x 1 matrix, where L is the number of levels
//*levels: Nx4 matrix with Bill's structure
//*theta: Covering radius shrinkage parameter
//*rootLevel: The integer index of the root level

//Outputs: *Edges: Mx3 list of edge indices; where M is the number of edges included
//The first two columns are the indices, and the third column is the relaxed distance
//*DeathTimes: Nx1 list of death times
void mexFunction(int nOutArray, mxArray *OutArray[], int nInArray, const mxArray *InArray[]) {  
	///////////////MEX INPUTS/////////////////
	if (nInArray < 2) {
		mexErrMsgTxt("Expecting point cloud as first input");
		return;
	}
	
	double* X;
	double* radii;
	int* levels;
	double theta;
	int rootLevel;
	
	const mwSize *dims;
	size_t N, D, NLevels;
	
	//X matrix
	int ndim = mxGetNumberOfDimensions(InArray[0]);
	if (ndim != 2) {
		mexErrMsgTxt("Expecting 2D matrix as first input");
		return;
	}
	dims = mxGetDimensions(InArray[0]);
	N = dims[0];
	D = dims[1];
	X = (double*)mxGetPr(InArray[0]);

	//Radii
	if (nInArray < 2) {
		mexErrMsgTxt("Expecting radii as second input");
		return;
	}
	dims = mxGetDimensions(InArray[1]);
	NLevels = dims[0];
	radii = (double*)mxGetPr(InArray[1]);
	
	//Levels
	if (nInArray < 3) {
		mexErrMsgTxt("Expecting levels as third input");
		return;
	}
	
	dims = mxGetDimensions(InArray[2]);
	if (dims[0] != N) {
		mexErrMsgTxt("Number of rows in levels != number of points");
		return;
	}
	if (dims[1] != 5) {
		mexErrMsgTxt("Should be 5 columns in levels");
		return;
	}
	levels = (int*)mxGetPr(InArray[2]);
	
	if (nInArray < 4) {
		mexErrMsgTxt("Expecting theta as fourth input");
		return;
	}
	theta = *((double*)mxGetPr(InArray[3]));
	
	if (nInArray < 5) {
		mexErrMsgTxt("Expecting root level as fifth input");
		return;
	}
	rootLevel = (int)(*((double*)mxGetPr(InArray[4])));
	
	///////////////ALGORITHM/////////////////
	
	vector<int> e1List;//Edges endpoint 1
	vector<int> e2List;//Edges endpoint 2
	vector<double> edList;//Edges relaxed distance
	
	//Calculate deletion times based on cover tree info
	double* ts = new double[N];
	for (int i = 0; i < N; i++) {
		int l = levels[i] - rootLevel;
		//mexPrintf("l = %i, t%i = ", l, i);
		ts[i] = 9*radii[l]/theta;//9*radius of parent
		//mexPrintf("%g\n", ts[i]);
	}
	
	//Check all paris of edges (slow version)
	for (int i = 0; i < N; i++) {
		for (int j = i+1; j < N; j++) {
			int e1, e2;
			double ed;
			getEdge(i, j, X, N, D, ts, &e1, &e2, &ed);
			if (e1 != -1 && e2 != -1) {
				e1List.push_back(e1);
				e2List.push_back(e2);
				edList.push_back(ed);
			}
		}
	}
	
	
	///////////////MEX OUTPUTS/////////////////
	size_t M = e1List.size() + N;
	double* edgeList = new double[M*3];
	//Add entries for points first
	for (size_t i = 0; i < N; i++) {
		edgeList[i] = i;
		edgeList[i+M] = i;
		edgeList[i+2*M] = 0;
	}
	//Now add edges between points
	for (size_t i = N; i < M; i++) {
		edgeList[i] = e1List[i-N];
		edgeList[i+M] = e2List[i-N];
		edgeList[i+2*M] = edList[i-N];
	}
	
	mwSize outdims[2];
	outdims[0] = M;
	outdims[1] = 3;
	OutArray[0] = mxCreateNumericArray(2, outdims, mxDOUBLE_CLASS, mxREAL);
	double* edgeListOut = (double*)mxGetPr(OutArray[0]);
	memcpy(edgeListOut, edgeList, M*3*sizeof(double));
	
	outdims[0] = N;
	outdims[1] = 1;
	OutArray[1] = mxCreateNumericArray(2, outdims, mxDOUBLE_CLASS, mxREAL);
	double* tsOut = (double*)mxGetPr(OutArray[1]);
	memcpy(tsOut, ts, N*sizeof(double));
	
	///////////////CLEANUP/////////////////
	delete[] ts;
	delete[] edgeList;
}
