addpath('BillsCode');
addpath('PCExamples');
NTHREADS=int32(4);

A.theta=.7;
A.numlevels=int32(20);
A.minlevel=int32(0);
A.NTHREADS=NTHREADS;
A.BLOCKSIZE=int32(32);

X = makeCircles();
X = X';
N = size(X, 2);
B = covertree(A, X);

%   int* levels;
%   int* parents;
%   int* numchildren;
%   int* childoffsets;
%   int* children;

rootLevel = B.outparams(2);

th = 0:0.1:2*pi+0.1;
plot(X(1, :), X(2, :), '.');
hold on;
numtext = {};
for kk = 1:N
   numtext{end+1} = sprintf('%i', B.levels(kk, 1) - rootLevel + 1);
end
text(X(1, :) + 0.001, X(2, :) + 0.001, numtext);


count = 0;
for ii = 1:N
    X1 = X(:, ii);
    NChildren = B.levels(ii, 3);
    if NChildren == 0
        continue;
    end
    for jj = (1:NChildren) + B.levels(ii, 4)
        child = B.levels(jj, 5)+1;
        X2 = X(:, child);
        plot([X1(1) X2(1)], [X1(2) X2(2)], 'r');
        count = count + 1;
    end
    
end
fprintf(1, '%i Links drawn\n', count);