init();%Initialize TDA tools

addpath('BillsCode');
addpath('PCExamples');
NTHREADS=int32(4);

A.theta=.5;
A.numlevels=int32(20);
A.minlevel=int32(0);
A.NTHREADS=NTHREADS;
A.BLOCKSIZE=int32(32);

%X = makeCircles(500);
X = rand(1000, 2);
N = size(X, 1);

%Step 1: Do the slow way, including all O(N^2) edges
tic;
I1 = rca1pc(X, 100);
time1 = toc;

%Step 2: Do the fast way, building a sparse edge list
tic;
%Build cover tree
B = covertree(A, X');
rootLevel = B.outparams(2);
B.levels(:, 1) = min(B.levels(:, 1), A.numlevels);
%Get sparse edge list from cover tree
[S, ts] = SlowSparseEdgeList(X, B.radii, B.levels, B.theta, rootLevel);
%Take care of 1-indexing
S(:, 1) = S(:, 1) + 1;
S(:, 2) = S(:, 2) + 1;
I2 = rca1mfscm(S, 100);
time2 = toc;



%Plot results
clf;
subplot(1, 2, 1);
hold on;
% for ii = 1:size(S, 1)
%     idx = [S(ii, 1), S(ii, 2)];
%     plot(X(idx, 1), X(idx, 2), 'b');
% end
plot(X(:, 1), X(:, 2), 'r.');

TotalEdges =  size(X, 1)*(size(X,1)-1)/2;
EdgesAdded = size(S, 1)-N;
title(sprintf('%i Verts: %i of %i Edges Added (%g %s)', size(X, 1), EdgesAdded, TotalEdges, 100*EdgesAdded/TotalEdges, '%'));


subplot(1, 2, 2);
scatter(I1(:, 1), I1(:, 2), 40, 'fill', 'r');
hold on;
scatter(I2(:, 1), I2(:, 2), 40, 'fill', 'b');
plot(I1(:, 1), I1(:, 2), 'rx');
legend({'Original', 'Approximate'});

plot([min(I1(:)) max(I1(:))], [min(I1(:)) max(I1(:))], 'r');
title(sprintf('%g sec vs %g sec', time1, time2));