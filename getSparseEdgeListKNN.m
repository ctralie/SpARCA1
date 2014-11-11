function [ S ] = getSparseEdgeListKNN( X, k )
    addpath('BillsCode');
    
    NTHREADS=int32(4);
    
    A.theta=.5;
    A.numlevels=int32(20);
    A.minlevel=int32(0);
    A.NTHREADS=NTHREADS;
    A.BLOCKSIZE=int32(32);

    N = size(X, 1);
    B = covertree(A, X');
    rootLevel = B.outparams(2);
    
    AllEdges = java.util.HashSet();
    
    for ii = 1:A.numlevels
        level = ii + rootLevel - 1;
        CentersIdx = 1:N;
        CentersIdx = CentersIdx(B.levels(:, 1) <= level);
        D = squareform(pdist(X(CentersIdx, :)));
        
        for jj = 1:size(D, 1)
            [~, idx] = sort(D(jj, :));
            knn = CentersIdx(idx(1:min(length(idx), k+1)));
            %Pick out all subsets of size 3 from the k nearest neighbors
            %and the vertex, inclusive
            for a = 1:length(knn)
                for b = a+1:length(knn)
                    for c = b+1:length(knn)
                        idx = sort([knn(a) knn(b) knn(c)]);
                        AllEdges.add(sprintf('%i,%i', idx(1), idx(2)));
                        AllEdges.add(sprintf('%i,%i', idx(2), idx(3)));
                        AllEdges.add(sprintf('%i,%i', idx(1), idx(3)));
                    end
                end
            end
        end
    end
    NEdges = AllEdges.size();
    AllEdges = AllEdges.toArray();
    S = zeros(NEdges + N, 3);
    S(1:N, 1:2) = repmat((1:N)', [1, 2]);
    for ii = 1:NEdges
        edge = str2num(AllEdges(ii));
        S(ii+N, 1) = edge(1);
        S(ii+N, 2) = edge(2);
        D = squareform(pdist(X(edge, :)));
        S(ii+N, 3) = D(2, 1);
    end
end