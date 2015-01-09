//Programmer: Chris Tralie
//Purpose: To extract a sparse edge list from a point cloud given a cover tree
//over that point cloud, to find cliques the slow way, and to compute dgm0 to dgmk
//using the naive reduction algorithm
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
#include <sstream>

using namespace std;

typedef struct sinf {
	string str;//TODO: Make this a pointer?
	size_t index;
	size_t k;//Number of vertices in the simplex
	double dist;
} SimplexInfo;

struct Simplex_DistComparator {
	bool operator()(SimplexInfo* s1, SimplexInfo* s2) const {
		if (s1->dist < s2->dist) {
			return true;
		}
		else if (s1->dist == s2->dist) {
			//Make sure to add the faces of a simplex before adding the simplex
			return s1->k < s2->k;
		}
		return false;
	}
};

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
void getEdgeRelaxedDist(int p, int q, double* X, int N, int D, double* ts, int* e1, int* e2, double* ed) {
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

//Return strings for all of the co-dimension 1 faces of the simplex 
//represented by str, and store them in the vector "cocliques"
void getCoDim1CliqueStrings(vector<string>& cocliques, string str) {
	int lasti = 0;
	for (size_t i = 0; i < str.size(); i++) {
		if (str[i] == '_') {
			//Omit the number between lasti and this index
			string str2 = "";
			if (lasti > 0) {
				str2 = str2 + str.substr(0, lasti);
				if (i + 1 < str.size()) {
					str2 = str2 + "_" + str.substr(i + 1);
				}
			}
			else {
				str2 = str.substr(i+1);
			}
			cocliques.push_back(str2);
			lasti = i;
		}
	}
}

//Check all subsets of a proximity list to find cliques by 
//using binary counting
void addCliquesFromProximityList(vector<int>& Ep, map<string, double>& EDists, map<string, SimplexInfo>* simplices, int MaxBetti) {
	if (Ep.size() == 0) {
		mexPrintf("Warning: Ep.size() = 0\n");
		return;
	}
	mexPrintf("Ep.size() = %i\n", Ep.size());
	mexEvalString("drawnow");
	stringstream ss;

	map<string, double>::iterator it;
	
	//Pre-copy pairwise edge lists from EDists
	//Cuts down on queries to the 
	int EpSize = (int)(Ep.size());
	double* dists = new double[EpSize*EpSize];
	int i1, i2, index;
	
	for (size_t i = 0; i < Ep.size(); i++) {
		for (size_t j = i+1; j < Ep.size(); j++) {
			ss.str("");
			ss << Ep[i] << "_" << Ep[j];
			it = EDists.find(ss.str());
			i1 = i*EpSize + j;
			i2 = j*EpSize + i;
			if (it == EDists.end()) {
				dists[i1] = -1;//TODO: Use something else for max distance?
				dists[i2] = -1;
			}
			else {
				dists[i1] = it->second;
				dists[i2] = it->second;
			}
		} 
	}
	
	//Now find all subsets up to size MaxBetti+2
	int MaxClique = min(MaxBetti+2, EpSize);
	vector<int> cliqueidx;
	while(true) {
		if (cliqueidx.size() < MaxClique) {
			if (cliqueidx.size() == 0) {
				cliqueidx.push_back(-1);//Base case
			}
			else {
				cliqueidx.push_back(cliqueidx[((int)cliqueidx.size()) - 1]);
			}
		}
		cliqueidx[((int)cliqueidx.size()) - 1]++;
		while (cliqueidx.size() > 1 && cliqueidx[((int)cliqueidx.size()) - 1] >= EpSize) {
			cliqueidx.pop_back();
			cliqueidx[((int)cliqueidx.size()) - 1]++;
		}
		if (cliqueidx.size() == 1 && cliqueidx[0] == EpSize) {
			break;
		}
		
		//Debugging output for mex to make sure I'm going through cliques properly
		/*ss.str("");
		ss << "(" << cliqueidx.size() << "): ";
		for (size_t i = 0; i < cliqueidx.size(); i++) {
			 ss << cliqueidx[i] << ", ";
		}
		ss << endl;
		mexPrintf("%s", ss.str().c_str());*/		
		
		//Now check every edge in the subset list to see
		//if it's a valid clique
		double maxDist = 0;
		bool validClique = true;
		for (size_t a = 0; a < cliqueidx.size(); a++) {
			for (size_t b = a + 1; b < cliqueidx.size(); b++) {
				i1 = cliqueidx[a];
				i2 = cliqueidx[b];
				index = i1*EpSize + i2;
				if (dists[index] == -1) {
					validClique = false;
				}
				maxDist = max(maxDist, dists[index]);
			}
			if (!validClique)
				break;
		}
		
		//If it's a valid clique, add it to the appropriate
		//sparse simplices map
		if (validClique) {
			vector<int> clique(cliqueidx.size());
			for (size_t k = 0; k < cliqueidx.size(); k++) {
				clique[k] = Ep[cliqueidx[k]];
			}
			//The indexing string for cliques is the index of all simplices
			//sorted in ascending order, separated by underscores
			sort(clique.begin(), clique.end());
			ss.str("");
			for (size_t k = 0; k < clique.size(); k++) {
				ss << clique[k];
				if ((int)k < ((int)clique.size()) - 1) {
					ss << "_";
				}
			}
			(simplices[((int)clique.size()) - 1])[ss.str()].index = 0;
			(simplices[((int)clique.size()) - 1])[ss.str()].dist = maxDist;
			(simplices[((int)clique.size()) - 1])[ss.str()].str = ss.str();
			(simplices[((int)clique.size()) - 1])[ss.str()].k = clique.size();
		}
	}
	delete[] dists;
}

//Inputs: *X: NxD point cloud matrix
//*radii: L x 1 matrix, where L is the number of levels
//*levels: Nx4 matrix with Bill's structure
//*theta: Covering radius shrinkage parameter
//*rootLevel: The integer index of the root level
//*MaxBetti: The maximum Betti numbers to compute (involving MaxK+2 simplices)

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
	int MaxBetti;
	
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
	
	if (nInArray < 6) {
		mexErrMsgTxt("Expecting MaxBetti as sixth input");
		return;
	}
	MaxBetti = (int)(*((double*)mxGetPr(InArray[5])));
	
	///////////////ALGORITHM/////////////////
	vector<int>* Ep = new vector<int>[N];//Proximity lists
	map<string, double> EDists;//Map from "v1_v2" to relaxed distance
	
	//Calculate deletion times based on cover tree info
	double* ts = new double[N];
	for (int i = 0; i < N; i++) {
		int l = levels[i] - rootLevel;
		//mexPrintf("l = %i, t%i = ", l, i);
		ts[i] = 9*radii[l]/theta;//9*radius of parent
		//mexPrintf("%g\n", ts[i]);
	}
	
	//Check all pairs of edges to find the sparse list (slow version)
	for (int i = 0; i < N; i++) {
		for (int j = i; j < N; j++) {
			int e1, e2;
			double ed;
			if (i == j) {
				ed = 0;//Include the vertex itself in the proximity list
				e1 = i;
				e2 = j;
			}
			else {
				getEdgeRelaxedDist(i, j, X, N, D, ts, &e1, &e2, &ed);
			}
			if (e1 != -1 && e2 != -1) {
				Ep[i].push_back(j);
				Ep[j].push_back(i);
				stringstream ss;
				ss << i << "_" << j;
				EDists[ss.str()] = ed;				
			}
		}
	}
	mexPrintf("Finished getting edge list\n");
	
	//Check all subsets of proximity lists to detect cliques
	map<string, SimplexInfo>* simplices = new map<string, SimplexInfo>[MaxBetti+2];//Index of simplices
	
	for (int i = 0; i < N; i++) {
		mexPrintf("%i of %i\n", i, N);
		mexEvalString("drawnow");
		addCliquesFromProximityList(Ep[i], EDists, simplices, MaxBetti);
	}
	
	//Sort the simplices in ascending order of distance
	for (int k = 0; k < MaxBetti+2; k++) {
		vector<SimplexInfo*> s;
		for (map<string,SimplexInfo>::iterator it = simplices[k].begin(); it != simplices[k].end(); it++) {
			s.push_back(&(it->second));
		}
		sort(s.begin(), s.end(), Simplex_DistComparator());
		for (size_t i = 0; i < s.size(); i++) {
			s[i]->index = i;
		}
	}
	
	//Debug print
	for (int k = 0; k < MaxBetti+2; k++) {
		mexPrintf("------------------------\nThere are %i %i-simplices\n", simplices[k].size(), k);
		for (map<string,SimplexInfo>::iterator it = simplices[k].begin(); it != simplices[k].end(); it++) {
			string str = it->first;
			mexPrintf("%s\n", str.c_str());
		}		
	}
	
	//Create the sparse boundary matrices for every level
	
	
	
	///////////////MEX OUTPUTS/////////////////
	
	/*mwSize outdims[2];
	outdims[0] = M;
	outdims[1] = 3;
	OutArray[0] = mxCreateNumericArray(2, outdims, mxDOUBLE_CLASS, mxREAL);
	double* edgeListOut = (double*)mxGetPr(OutArray[0]);
	memcpy(edgeListOut, edgeList, M*3*sizeof(double));
	
	outdims[0] = N;
	outdims[1] = 1;
	OutArray[1] = mxCreateNumericArray(2, outdims, mxDOUBLE_CLASS, mxREAL);
	double* tsOut = (double*)mxGetPr(OutArray[1]);
	memcpy(tsOut, ts, N*sizeof(double));*/
	
	///////////////CLEANUP/////////////////
	delete[] ts;
	delete[] Ep;
	delete[] simplices;
}
