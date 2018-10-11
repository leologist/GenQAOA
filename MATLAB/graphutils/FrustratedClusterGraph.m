function [edges, AdjMat] = FrustratedClusterGraph(n, a, b)
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

AdjMat = zeros(4*n);

for ind = 1:n
    for ind2 = 1:n
        AdjMat(ind, n+ind2) = 1;
        AdjMat(2*n+ind, 3*n+ind2) = 1;
    end
end

for ind = 1:a
    AdjMat(ind, 2*n+ind) = 1;
end

for ind = 1:b
    AdjMat(ind, 3*n+ind) = 1;
end

[row, col] = find(AdjMat);
edges = [row, col];

AdjMat = AdjMat + AdjMat';



end

