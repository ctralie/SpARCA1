MEX = mex #Using matlab
#MEX = mkoctfile --mex#Using octave
#Change the path below to match your matlab path
MEXINCLUDE = -I/usr/local/MATLAB/R2014b/extern/include/


#LIBS = -lcudart -lcublas

all: NaiveSparseCliqueReduction SlowSparseEdgeList GreedyKCenter NaiveGreedyKCenter

NaiveSparseCliqueReduction: NaiveSparseCliqueReduction.cpp
	$(MEX) -g NaiveSparseCliqueReduction.cpp $(MEXINCLUDE)

SparseEdgeList: SparseEdgeList.cpp
	$(MEX) -g SparseEdgeList.cpp $(MEXINCLUDE)

SlowSparseEdgeList: SlowSparseEdgeList.cpp
	$(MEX) -g SlowSparseEdgeList.cpp $(MEXINCLUDE)

GreedyKCenter: GreedyKCenter.cpp
	$(MEX) -g GreedyKCenter.cpp $(MEXINCLUDE)

NaiveGreedyKCenter: NaiveGreedyKCenter.cpp
	$(MEX) -g NaiveGreedyKCenter.cpp

clean:
	rm -f *.mexa64
