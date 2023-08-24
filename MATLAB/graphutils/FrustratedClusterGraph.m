function [edges, AdjMat] = FrustratedClusterGraph(n, a, b, randomize)
%FrustratedClusterGraph generate a frustrated cluster graph that should be
%bad for simulated annealing approach to find MAXCUT
%   [edges, AdjMat] = FrustratedClusterGraph(n, a, b)
%   generates a graph with 4*n nodes with the following structure
%       - nodes are divided into four clusters of equal sizes: A, B, C, D
%       - complete bipartite graph between A and B
%       - complete bipartite graph between C and D
%       ~ there are a edges between A and C
%       - there are b edges between A and D
%
%   This is bad for simulated annealing because there are two local maxmima
%   of MAXCUT objective function separated by large Hamming distance.

if nargin <= 3
    randomize = false;
end

AdjMat = zeros(4*n);

for ind = 1:n
    for ind2 = 1:n
        AdjMat(ind, n+ind2) = 1;
        AdjMat(2*n+ind, 3*n+ind2) = 1;
    end
end

% pairings = nchoosek(1:n,2);
[p1,p2] = meshgrid(1:n,1:n);
pairings = [p1(:), p2(:)];
if randomize
    pairings = pairings(randperm(size(pairings,1)),:);
end

for ind = 1:a
    AdjMat(pairings(ind,1), 2*n+pairings(ind,2)) = 1;
end

if randomize
    pairings = pairings(randperm(size(pairings,1)),:);
end

for ind = 1:b
%     AdjMat(ind, 3*n+ind) = 1;
    AdjMat(pairings(ind,1), 3*n+pairings(ind,2)) = 1;
end

[row, col] = find(AdjMat);
edges = [row, col];

AdjMat = AdjMat + AdjMat';



end

