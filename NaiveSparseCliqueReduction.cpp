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

#define VERBOSE 0

typedef struct sinf {
	string str;//TODO: Make this a pointer?
	size_t index;
	size_t k;//Number of vertices in the simplex
	double dist;
} SimplexInfo;

typedef struct bdt {
	double birth;
	double death;
} BDTimes;

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

//Print a sparse matrix (for debugging)
void printMatrix(vector<int>* M, int n, int m) {
	int* colPointers = new int[m];
	for (int i = 0; i < m; i++) {
		colPointers[i] = 0;
	}
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < m; j++) {
			int elem = 0;
			if (colPointers[j] < M[j].size()) {
				if (i == M[j][colPointers[j]]) {
					elem = 1;
					colPointers[j]++;
				}
			}
			mexPrintf("%i ", elem);
		}
		mexPrintf("\n");
	}
	delete[] colPointers;
}

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

//Function returns the original (unwarped) metric
void getEdgeUnwarpedDist(int p, int q, double* X, int N, int D, double* ts, int* e1, int* e2, double* ed) {
	double dpq = sqrt(getSqrDist(X, N, D, p, q));
	*e1 = p;
	*e2 = q;
	*ed = dpq;
}

//Return strings for all of the co-dimension 1 faces of the simplex 
//represented by str, and store them in the vector "cocliques"
void getCoDim1CliqueStrings(vector<string>& cocliques, string str) {
	size_t lasti = 0;
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
	//Add the last string resulting from omitting the last element
	cocliques.push_back(str.substr(0, lasti));
}

//Check all subsets of a proximity list to find cliques by 
//using binary counting
void addCliquesFromProximityList(vector<int>& Ep, map<string, double>& EDists, map<string, SimplexInfo>* simplices, int MaxBetti) {
	if (Ep.size() == 0) {
		mexPrintf("Warning: Ep.size() = 0\n");
		return;
	}
	if (VERBOSE) {
		mexPrintf("Ep.size() = %i\n", Ep.size());
		mexEvalString("drawnow");
	}
	stringstream ss;

	map<string, double>::iterator it;
	
	/***Pre-copy pairwise edge lists from EDists***/
	//Cuts down on queries to the distance map
	int EpSize = (int)(Ep.size());
	double* dists = new double[EpSize*EpSize];
	size_t i1, i2, index;
	
	for (size_t i = 0; i < Ep.size(); i++) {
		if (VERBOSE) {
			mexPrintf("%i ", Ep[i]);
		}
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
	if (VERBOSE) {
		mexPrintf("\n");
		mexEvalString("drawnow");
	}
	
	/***Now find all subsets up to size MaxBetti+2***/
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

//Do linear time addition of two sorted columns in the sparse matrix
//This method assumes that col1 and col2 are in sorted order by TDASimplex.id (the filtration)
void addColToColMod2(vector<int>& col1, vector<int>& col2) {
	int i1 = 0, i2 = 0;
	vector<int> out;
	while (i1 < (int)col1.size() && i2 < (int)col2.size()) {
		int id1 = col1[i1];
		int id2 = col2[i2];
		if (id1 == id2) {
			i1++;
			i2++;
			//Do nothing; they cancel out mod2 if they are the same
		}
		else if (id1 < id2) {
			//Add the element from col1 first and move down one element on col1
			out.push_back(col1[i1]);
			i1++;
		}
		else if (id2 < id1) {
			//Add the element from col2 first and move down one element on col2
			out.push_back(col2[i2]);
			i2++;
		}
		if (i1 == (int)col1.size()) {
			//Add the rest of the elements from column 2 if I'm through column 1
			while (i2 < (int)col2.size()) {
				out.push_back(col2[i2]);
				i2++;
			}
		}
		if (i2 == (int)col2.size()) {
			//Add the rest of the elements from column 1 if I'm through column 2
			while (i1 < (int)col1.size()) {
				out.push_back(col1[i1]);
				i1++;
			}
		}
	}
	//Overwrite col2 with the result
	col2.clear();
	col2.insert(col2.begin(), out.begin(), out.end());
}

//Add the column "col" to every subsequent column in M that contains the low 
//element of col in the matrix M, which has m columns
//This method assumes that the columns M are in sorted order
void addLowElementToOthers(vector<int>* M, int col, int m) {
	assert(col >= 0 && col < m);
	assert(M[col].size() > 0);
	int low = (M[col])[M[col].size()-1];//The low element is the last element in sorted order
	for (int j = col+1; j < m; j++) {
		//Check to see if this column has the low element
		for (size_t i = 0; i < M[j].size(); i++) {
			if (M[j][i] == low) {
				//This column does contain the low element so add M[col] to it
				addColToColMod2(M[col], M[j]);
				break;
			}
		}
	}
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
	mexPrintf("\n\n=== Running Chris Tralie's sparse clique reduction for higher homology ====\n");
	vector<int>* Ep = new vector<int>[N];//Proximity lists
	map<string, double> EDists;//Map from "v1_v2" to relaxed distance
	
	/***Calculate deletion times based on cover tree info***/
	double* ts = new double[N];
	for (int i = 0; i < N; i++) {
		int l = levels[i] - rootLevel;
		ts[i] = 9*radii[l]/theta;//9*radius of parent
	}
	
	/***Check all pairs of edges to find the sparse list (slow version)***/
	for (int i = 0; i < N; i++) {
		for (int j = i; j < N; j++) {
			int e1, e2;
			double ed;
			if (i == j) {
				//Include the vertex itself in the proximity list
				Ep[i].push_back(i);
				ed = 0;
				stringstream ss;
				ss << i << "_" << j;
				EDists[ss.str()] = ed;			
			}
			else {
				getEdgeRelaxedDist(i, j, X, (int)N, (int)D, ts, &e1, &e2, &ed);
				//getEdgeUnwarpedDist(i, j, X, (int)N, (int)D, ts, &e1, &e2, &ed);
				if (e1 != -1 && e2 != -1) {					
					Ep[i].push_back(j);
					Ep[j].push_back(i);
					stringstream ss;
					ss << i << "_" << j;
					EDists[ss.str()] = ed;				
				}
			}

		}
	}
	mexPrintf("Finished getting edge list\n");
	mexPrintf("Checking proximity list subsets to find cliques...");
	/***Check all subsets of proximity lists to detect cliques***/
	map<string, SimplexInfo>* simplices = new map<string, SimplexInfo>[MaxBetti+2];//Index of simplices
	for (int i = 0; i < N; i++) {
		if (i % 50 == 0)
			mexPrintf("\n");
		mexPrintf(".", i, N);
		mexEvalString("drawnow");
		addCliquesFromProximityList(Ep[i], EDists, simplices, MaxBetti);
	}
	mexPrintf("\n\n");
	
	/***Sort the simplices in ascending order of distance***/
	SimplexInfo*** simplexInfoSorted = new SimplexInfo**[MaxBetti+2];
	for (int k = 0; k < MaxBetti+2; k++) {
		vector<SimplexInfo*> s;
		for (map<string,SimplexInfo>::iterator it = simplices[k].begin(); it != simplices[k].end(); it++) {
			s.push_back(&(it->second));
		}
		sort(s.begin(), s.end(), Simplex_DistComparator());
		simplexInfoSorted[k] = new SimplexInfo*[s.size()];
		for (size_t i = 0; i < s.size(); i++) {
			s[i]->index = i;
			simplexInfoSorted[k][i] = s[i];
		}
	}
	
	/***Print number of simplices of each order and report savings***/
	for (int k = 0; k < MaxBetti+2; k++) {
		//Compute number of combinations explicitly
		int N2 = 0;
		int numer = 1;
		for (int i = (int)N; i > (int)N - k - 1; i--) {
			numer = numer*i;
		}
		int denom = 1;
		for (int i = 1; i <= k+1; i++) {
			denom = denom*i;
		}
		N2 = numer/denom;
	
		mexPrintf("------------------------\nThere are %i %i-simplices (%i max possible: %g%%)\n", simplices[k].size(), k, N2, 100*((double)simplices[k].size())/((double)N2));
	}

	
	/***Create the sparse boundary matrices for every level***/
	vector<int>** Ms = new vector<int>*[MaxBetti+2];
	for (int k = 0; k < MaxBetti+2; k++) {
		Ms[k] = new vector<int>[simplices[k].size()];
		if (k == 0) {
			continue; //The zero boundary matrix is a dummy one
		}
		for (int j = 0; j < simplices[k].size(); j++) {
			vector<string> cocliques;
			getCoDim1CliqueStrings(cocliques, simplexInfoSorted[k][j]->str);
			for (size_t a = 0; a < cocliques.size(); a++) {
				int cofaceIndex = simplices[k-1][cocliques[a]].index;
				Ms[k][j].push_back(cofaceIndex);
			}
			//Keep the invariant that sparse columns elements are in order
			sort(Ms[k][j].begin(), Ms[k][j].end());
		}	
	}
	
	/***Reduce the matrices to compute the persistence diagrams***/
	//Index the homology classes by the index of their birthing simplex
	map<int, BDTimes>* Is = new map<int, BDTimes>[MaxBetti+1];
	//Assume for now all points are born at 0
	for (int i = 0; i < N; i++) {
		Is[0][i].birth = 0;
		Is[0][i].death = -1;
	}
	//Now reduce the rest of the boundary matrices
	mexPrintf("\n\n");
	for (int k = 1; k < MaxBetti+2; k++) {
		mexPrintf("Reducing %i-dimensional boundary matrix...\n", k);
		mexEvalString("drawnow");
		for (size_t j = 0; j < simplices[k].size(); j++) {
			//mexPrintf(" %g",  simplexInfoSorted[k][j]->dist);
			if (Ms[k][j].size() > 0) {
				//mexPrintf("(1)");
				addLowElementToOthers(Ms[k], (int)j, simplices[k].size());
				//Rank of Mk increases by one, so a (k-1)-class is killed
				//Pair with the (k-1) co-face that most recently participated in a birth
				//(the element at the bottom of the column)
				int killIndex = Ms[k][j][Ms[k][j].size()-1];
				map<int, BDTimes>::iterator iter = Is[k-1].find(killIndex);
				if (iter == Is[k-1].end()) {
					mexPrintf("Error: Trying to pair a death with a simplex that never gave rise to a birth\n");
					continue;
				}
				BDTimes* I = &(iter->second);
				if (I->death != -1) {
					mexPrintf("Warning: Trying to kill class that's already been killed at %g\n", I->death);
				}
				I->death = simplexInfoSorted[k][j]->dist;
				if (VERBOSE) {
					mexPrintf("%iD class born with simplex %s ", k-1, simplexInfoSorted[k-1][killIndex]->str.c_str());
					mexPrintf("killed by simplex %s ", simplexInfoSorted[k][j]->str.c_str());
					mexPrintf("(dist %g)\n", I->death);
				}
				if (I->death == I->birth) {
					//Don't include classes that are born and die instantly
					Is[k-1].erase(killIndex);
				}
			}
			else {
				//Kernel of Mk increases by one, so a k-class is born
				//(Ignore births for the last class)
				if (k <= MaxBetti) {
					BDTimes bd;
					bd.birth = simplexInfoSorted[k][j]->dist;
					bd.death = -1;
					(Is[k])[j] = bd;
				}
			}
		}
	}	
	
	///////////////MEX OUTPUTS/////////////////
	/***Output a cell array holding all of the persistence diagrams***/
	mwSize outDims[2];
	outDims[0] = MaxBetti+1;
	outDims[1] = 1;
	OutArray[0] = mxCreateCellArray(2, outDims);
	mxArray* cellArrayPtr = OutArray[0];
	for (int k = 0; k < MaxBetti+1; k++) {
		mxArray* I = mxCreateDoubleMatrix(Is[k].size(), 2, mxREAL);
		double* IArray = mxGetPr(I);
		int i = 0;
		for (map<int,BDTimes>::iterator it = Is[k].begin(); it != Is[k].end(); it++) {
			IArray[i] = it->second.birth;
			IArray[i+Is[k].size()] = it->second.death;
			i++;
		}
		//Place the matrix into the cell array
		mxSetCell(cellArrayPtr, k, I);
	}

	
	///////////////CLEANUP/////////////////
	delete[] Is;
	for (int i = 0; i < MaxBetti+2; i++) {
		delete[] Ms[i];
	}
	delete[] Ms;
	for (int i = 0; i < MaxBetti+2; i++) {
		delete[] simplexInfoSorted[i];
	}
	delete[] simplexInfoSorted;
	delete[] simplices;	
	delete[] Ep;
	delete[] ts;
}
