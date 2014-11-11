addpath('BillsCode');
addpath('PCExamples');
NTHREADS=int32(4);

A.theta=.5;
A.numlevels=int32(20);
A.minlevel=int32(0);
A.NTHREADS=NTHREADS;
A.BLOCKSIZE=int32(32);

X = makeCircles();
X = X';
N = size(X, 2);
B = covertree(A, X);

rootLevel = B.outparams(2);

colors = rand(N, 3);

th = 0:0.1:2*pi+0.1;
for ii = 1:10
    clf;
    level = ii + rootLevel - 1;
    RPack = B.radii(ii);
    RCover = (1/(1-A.theta))*RPack;%Make sure radius covers
    plot(X(1, :), X(2, :), 'r.');
    hold on;

    centers = getSubtreesAtLevel(B, level);
    for jj = 1:length(centers)
        idx = centers{jj}(1);
        plot(X(1, idx) + RPack*cos(th), X(2, idx) + RPack*sin(th), 'b');
        plot(X(1, idx) + RCover*cos(th), X(2, idx) + RCover*sin(th), 'r');
        scatter(X(1, idx), X(2, idx), 20, colors(jj, :));
        plot(X(1, centers{jj}), X(2, centers{jj}), '.', 'color', colors(jj, :));
    end
    
    xlim([min(X(1, :)) - B.radii(3), max(X(1, :)) + B.radii(3)]);
    ylim([min(X(2, :)) - B.radii(3), max(X(2, :)) + B.radii(3)]);
    axis equal;
    print('-dpng', '-r300', sprintf('%i.png', ii));
end