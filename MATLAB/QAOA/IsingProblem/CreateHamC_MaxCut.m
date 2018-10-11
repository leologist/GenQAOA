function HamC = CreateHamC_MaxCut(N, wGraph)
%   HamC = CreateHamC(N, wGraph)
%
%   This function creates the Hamiltonian for MaxCut objective
%   function, HamC = \sum_{ij} w_{ij} (I - Z_i Z_j)/2
%   We only output the diagonal (HamC is 2^N x 1 vector)
%
%   wGraph is in the format of (i,j,wij) as columns. When 
%   the third column is missing, we assume wij = 1

numEdges = size(wGraph,1);
if size(wGraph,2) < 3
    wGraph = [wGraph, ones(numEdges, 1)];
end
    
HamC = zeros(2^N,1);
for ind = 1:numEdges
    HamC = HamC + wGraph(ind,3)*(1-kronSz(N,wGraph(ind,1), wGraph(ind,2)))/2; % 1/2*(1-sz_i sz_j)
end


function [out] = kronSz(N,i,j)
%kronSz makes a kronecker product of N spin-1/2 spin matrices, representing  Sz^i Sz^j.
% Only the diagonal part are outputed since they are Sz matrices. 
% N: total number of sites. i,j are the sites where spin matrices are Sz. 
% All other sites have identity spin matrix. 

ind = i;
if i>=j, i = j; j = ind; end % swap i and j is i>=j
    
sz = [1; -1];
IdL=ones(2^(i-1),1); % i-1 is the number of sites on the left
IdM=ones(2^(j-i-1),1); % in the middle
IdR=ones(2^(N-j),1); % on the right

outT1= kron(IdL, sz);
outT2= kron(outT1, IdM);
outT3= kron(outT2, sz);
out= kron(outT3, IdR);