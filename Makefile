MEX = mex
MEXINCLUDE = -I/usr/local/MATLAB/R2013b/extern/include/


#LIBS = -lcudart -lcublas

all: GreedyKCenter NaiveGreedyKCenter

GreedyKCenter: GreedyKCenter.cpp
	$(MEX) -g GreedyKCenter.cpp $(MEXINCLUDE)

NaiveGreedyKCenter: NaiveGreedyKCenter.cpp
	$(MEX) -g NaiveGreedyKCenter.cpp

clean:
	rm -f *.mexa64
