function HamC = CreateHamC_Ising(N, h_i, J_ij)
%CreateHamC_Ising
%
%   HamC = \sum_i h_i Z_i + \sum_{<ij>} J_ij Z_i Z_j
%   Usage:
%       HamC = CreateHamC_Ising(N, h_i, J_ij)

if isempty(h_i)
    h_i = zeros(N, 1);
elseif numel(h_i) ~= N
    error('Invalid input h')
end
if size(J_ij, 2) == 2
    J_ij = [J_ij, ones(size(J_ij,1),1)];
end

HamC = krondist([1;-1], N, h_i);

for ind = 1:size(J_ij, 1)
    HamC = HamC + J_ij(ind, 3) * kronSz(N, J_ij(ind,1), J_ij(ind, 2));
end


end


function out = kronSz(N,i,j)
%kronSz makes a kronecker product of N spin-1/2 spin matrices, representing  Sz^i Sz^j.
% Only the diagonal part are outputed since they are Sz matrices. 
% N: total number of sites. i,j are the sites where spin matrices are Sz. 
% All other sites have identity spin matrix. 

ind = i;
if i>=j, i = j; j = ind; end % swap i and j is i>=j
    
if i == j
    out = ones(2^N,1);
    return
end
sz = [1; -1];
IdL=ones(2^(i-1),1); % i-1 is the number of sites on the left
IdM=ones(2^(j-i-1),1); % in the middle
IdR=ones(2^(N-j),1); % on the right

outT1= kron(IdL, sz);
outT2= kron(outT1, IdM);
outT3= kron(outT2, sz);
out= kron(outT3, IdR);

end