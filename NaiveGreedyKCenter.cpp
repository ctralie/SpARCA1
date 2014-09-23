//Programmer: Chris Tralie
//Purpose: To implement the naive O(n^2) algorithm to perform the
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
#include <map>


using namespace std;

//Replace x with min(x, y)
void replaceWithMin(double* x, double*y, int N) {
	for (int i = 0; i < N; i++) {
		if (y[i] < x[i]) {
			x[i] = y[i];
		}
	}
} 

//Return max and argmax of the array (x) by reference
void max(double* x, int N, double* max, int* index) {
	if (N <= 0) {
		*max = 0;
		*index = 0;
	}
	for (int i = 0; i < N; i++) {
		if (x[i] > *max) {
			*max = x[i];
			*index = i;
		}
	}
}

//Inputs: D (NxN distance matrix)
//Outputs: centers (1xN "permutation" of centers)
//rads (1xN matrix of radii needed to cover points at each level)
void mexFunction(int nOutArray, mxArray *OutArray[], int nInArray, const mxArray *InArray[]) {  
	///////////////MEX INPUTS/////////////////
	if (nInArray < 1) {
		mexErrMsgTxt("Distance Matrix Input Required");
		return;
	}
	
	const mwSize *dims;
	size_t rowsD, colsD;
	int ndim = mxGetNumberOfDimensions(InArray[0]);
	if (ndim != 2) {
		mexErrMsgTxt("Expecting 2D matrix as first input");
		return;
	}
	dims = mxGetDimensions(InArray[0]);
	rowsD = dims[0];
	colsD = dims[1];
	if (rowsD != colsD) {
		mexErrMsgTxt("Number of rows not equal to number of columns in distance matrix");
		return;
	}
	
	///////////////ALGORITHM/////////////////
	double* D = (double*)mxGetPr(InArray[0]);
	size_t N = rowsD;
	
	double* centers = new double[N];//Index of the centers chosen in the order they are chosen
	double* rads = new double[N];//Radius of the balls needed around the (i <= k) center needed
	//to cover all points at level k
	double* dists = new double[N];//Closest distances of all points to the centers chosen so far
	
	centers[0] = 1;//Make the first center the first point in the list
	//to keep things deterministic
	memcpy(dists, D, sizeof(double)*N);//Therefore the distances to the centers start
	//of as being the distances to the first point
	for (size_t i = 1; i < N; i++) {
		int index = 0;
		double R = 0;
		max(dists, N, &R, &index);
		centers[i] = (double)(index+1);//Matlab is 1-indexed
		rads[i-1] = R;
		//A new center has been added so update the distance of points that
		//may now be closer to this center than they were to previous centers
		replaceWithMin(dists, D + N*index, N);
	}
	rads[N - 1] = 0;//Doesn't matter what the radius at the last level is because it
	//includes all points
	
	///////////////MEX OUTPUTS/////////////////
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
	delete[] rads;
	delete[] dists;
}
