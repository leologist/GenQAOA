function [edges, AdjMat] = HammingGraph(d,h)
%This function generates the binary Hamming graph with 2^d vertices, and Hamming distance h
% see the book "Approximation Algorithms and Semidefinite Programming", page 149 for more details

N = 2^d; % number of vertices 
vertex = dec2bin(0:N-1,d);
vertex = double(vertex); % note: double converts '0' to 48 and '1' to 49. 
HamDist = round(d*pdist(vertex, 'Hamming'));
HamDist = squareform(HamDist); 
AdjMat = (HamDist==h);
edges = triu(AdjMat); % use only the upper triangular part 
[row, col] = find(edges);
edges = [row, col, ones(numel(row),1)]; % the last columns are weights = 1 
