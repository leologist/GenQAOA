function HamC = CreateHamNNInt(N, wGraph)
%CreateHamNNInt generate diagonal Hamiltonian H = \sum_{ij} w_ij n_i n_j


numEdges = size(wGraph,1);
if size(wGraph,2) < 3
    wGraph = [wGraph, ones(numEdges, 1)];
end

HamC = 0;

for ind = 1:numEdges
    HamC = HamC + wGraph(ind,3)*kronNz(N, wGraph(ind,1), wGraph(ind,2));
end


end

function [out] = kronNz(N,i,j)
%kronSz makes a kronecker product of N spin-1/2 spin matrices, representing  Sz^i Sz^j.
% Only the diagonal part are outputed since they are Sz matrices. 
% N: total number of sites. i,j are the sites where spin matrices are Sz. 
% All other sites have identity spin matrix. 

if i > j
    out = kronNz(N, j, i);
    return
end

nz = [0; 1];
IdL=ones(2^(i-1),1); % i-1 is the number of sites on the left
IdM=ones(2^(j-i-1),1); % in the middle
IdR=ones(2^(N-j),1); % on the right

outT1= kron(IdL, nz);
outT2= kron(outT1, IdM);
outT3= kron(outT2, nz);
out= kron(outT3, IdR);

end