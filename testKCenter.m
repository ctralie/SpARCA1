N = 10;
X = randn(N, 2);
D = squareform(pdist(X));
[centers, rads] = NaiveGreedyKCenter(D);

theta = 0:0.1:2*pi+0.1;
buffer = mean(D(:));
limX = [min(X(:, 1)) - buffer, max(X(:, 1)) + buffer];
limY = [min(X(:, 2)) - buffer, max(X(:, 2)) + buffer];
for k = 1:N;
    clf;
    plot(X(:, 1), X(:, 2), 'r.');
    hold on;
    R = rads(k);
    for ii = 1:k
       C = centers(ii);
       plot(X(C, 1) + R*cos(theta), X(C, 2) + R*sin(theta));
       scatter(X(C, 1), X(C, 2));
    end
    scatter(X(C, 1), X(C, 2), 100, 'k', 'fill');
    xlim(limX);
    ylim(limY);
    print('-dpng', '-r100', sprintf('%i.png', k));
end