% %%  Testing 2 circles case to make sure DGM1 is correct
% init();%Initialize TDA tools
% addpath('BillsCode');
% addpath('PCExamples');
% 
% NTHREADS=int32(4);
% A.theta=.5;
% A.numlevels=int32(20);
% A.minlevel=int32(0);
% A.NTHREADS=NTHREADS;
% A.BLOCKSIZE=int32(32);
% 
% n = 100;
% X = makeCircles(n);
% N = size(X, 1);
% 
% %Build cover tree
% B = covertree(A, X');
% rootLevel = B.outparams(2);
% B.levels(:, 1) = min(B.levels(:, 1), A.numlevels);
% 
% % % Do RCA1 on the sparse edge list
% % % Get sparse edge list from cover tree
% % [S, ts] = SlowSparseEdgeList(X, B.radii, B.levels, B.theta, rootLevel);
% % %Take care of 1-indexing
% % S(:, 1) = S(:, 1) + 1;
% % S(:, 2) = S(:, 2) + 1;
% % [I2Orig, I20Orig] = rca1mfscm(S, 100);
% % IsOrig = {I20Orig, I2Orig};
% tic
% [I2OrigFull, I20OrigFull] = rca1pc(X, 1e9);
% IsOrigFull = {I20OrigFull, I2OrigFull};
% toc
% 
% save(sprintf('PCExamples/Circles%i.mat', n), 'A', 'X', 'B', 'rootLevel');
% 
% % load('PCExamples/Circles25.mat');
% tic
% Is = NaiveSparseCliqueReduction(X, B.radii, B.levels, B.theta, rootLevel, 1);
% toc
% 
% clf;
% for ii = 1:length(Is)
%     subplot(1, 2, ii);
%     plotpersistencediagram(Is{ii});
%     if (ii <= 2) 
%         hold on;
%         scatter(Is{ii}(:, 1), Is{ii}(:, 2), 30, 'b', 'fill');
%         plot(IsOrigFull{ii}(:, 1), IsOrigFull{ii}(:, 2), 'rx');%Ground truth from RCA1
%     end
%     title(sprintf('DGM %i (%i Classes)', ii-1, size(Is{ii}, 1)));
% end

%% Testing 2-torus case to make sure up to 2D homology is correct
init();%Initialize TDA tools
addpath('BillsCode');
addpath('PCExamples');

R1 = 4;
R2 = 1.5;
N = 200;
X = make2Torus(R1, R2, N);
%X = make2Sphere(1, N);
N = size(X, 1);

NTHREADS=int32(4);
A.theta=.1;
A.numlevels=int32(N);
A.minlevel=int32(0);
A.NTHREADS=NTHREADS;
A.BLOCKSIZE=int32(20);

%Build cover tree
B = covertree(A, X');
rootLevel = B.outparams(2);
B.levels(:, 1) = min(B.levels(:, 1), A.numlevels);

save(sprintf('PCExamples/2Torus%i.mat', N), 'A', 'X', 'B', 'rootLevel');
%load('PCExamples/2Sphere121.mat');


%Run up to DGM1 on the full unwarped rips filtration
tic
[I2OrigFull, I20OrigFull] = rca1pc(X, 1e9);
IsOrigFull = {I20OrigFull, I2OrigFull};
toc

tic
Is = NaiveSparseCliqueReduction(X, B.radii, B.levels, B.theta, rootLevel, 2);
toc

clf;
subplot(2, 2, 1);
plot3(X(:, 1), X(:, 2), X(:, 3), '.');
axis equal;
title(sprintf('2 Torus R1 = %g, R2 = %g, %i Points', R1, R2, N));
%title(sprintf('2 Sphere %i points', N));
for ii = 1:length(Is)
    subplot(2, 2, ii+1);
    if isempty(Is{ii})
        title(sprintf('DGM %i (Empty)', ii-1));
        continue;
    end
    plotpersistencediagram(Is{ii});
    hold on;
    scatter(Is{ii}(:, 1), Is{ii}(:, 2), 30, 'b', 'fill');
    if ii <= 2
        plot(IsOrigFull{ii}(:, 1), IsOrigFull{ii}(:, 2), 'rx');%Ground truth from RCA1
    end
    title(sprintf('DGM %i (%i Classes)', ii-1, size(Is{ii}, 1)));
end

save('Torus.mat', 'Is', 'IsOrigFull');

setpref('Internet', 'E_mail', 'labmailciemas3418@gmail.com');
setpref('Internet', 'SMTP_Server', 'smtp.gmail.com');
setpref('Internet', 'SMTP_Username', 'labmailciemas3418@gmail.com');
setpref('Internet', 'SMTP_Password', 'r3yn0ldsL@b');
props = java.lang.System.getProperties;
props.setProperty('mail.smtp.auth','true');
props.setProperty('mail.smtp.socketFactory.class', 'javax.net.ssl.SSLSocketFactory');
props.setProperty('mail.smtp.socketFactory.port','465');
sendmail('chris.tralie@gmail.com', sprintf('Finished 2 Torus R1 = %g, R2 = %g, %i Points', R1, R2, N));
