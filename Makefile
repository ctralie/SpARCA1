MEX = mex
MEXINCLUDE = -I/usr/local/MATLAB/R2013b/extern/include/


#LIBS = -lcudart -lcublas

all: GreedyKCenter NaiveGreedyKCenter

GreedyKCenter: GreedyKCenter.cpp
	$(MEX) -g GreedyKCenter.cpp $(LIBS)

NaiveGreedyKCenter: NaiveGreedyKCenter.cpp
	$(MEX) -g NaiveGreedyKCenter.cpp $(LIBS)

clean:
	rm -f *.mexa64
