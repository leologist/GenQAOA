function [edges, AdjMat] = ErdosRenyiGraph(N, p)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

AdjMat = rand(N,N) <= p;
AdjMat = triu(AdjMat, 1);

[r,c] = find(AdjMat);
edges = [r, c];

AdjMat = AdjMat + AdjMat';

end

