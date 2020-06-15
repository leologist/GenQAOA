function [edges, AdjMat] = ErdosRenyiGraph(N, p)
%ErdosRenyiGraph generate an Erdos-Renyi graph of N vertices and edge
%   probability p
%

AdjMat = rand(N,N) <= p;
AdjMat = triu(AdjMat, 1);

[r,c] = find(AdjMat);
edges = [r, c];

AdjMat = AdjMat + AdjMat';

end

