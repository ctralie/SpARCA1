init();
addpath('PCExamples');
X = rand(400, 2);

tic;
S = getSparseEdgeListKNN(X, 4);
time2 = toc;

clf;
subplot(1, 2, 1);
plot(X(:, 1), X(:, 2), 'r.');
hold on;
for ii = 1:size(S, 1)
    idx = [S(ii, 1), S(ii, 2)];
    plot(X(idx, 1), X(idx, 2), 'b');
end
title(sprintf('%i Verts: %i of %i Edges Added', size(X, 1), size(S, 1) - size(X, 1), size(X, 1)*(size(X,1)-1)/2));

tic;
I1 = rca1pc(X, 100);
time1 = toc;
tic;
I2 = rca1mfscm(S, 100);
time2 = time2 + toc;

subplot(1, 2, 2);
scatter(I1(:, 1), I1(:, 2), 40, 'fill', 'r');
hold on;
scatter(I2(:, 1), I2(:, 2), 40, 'fill', 'b');
legend({'Original', 'Approximate'});

plot([min(I1(:)) max(I1(:))], [min(I1(:)) max(I1(:))], 'r');
title(sprintf('%g sec vs %g sec', time1, time2));