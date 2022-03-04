function [HCsym, HBsym, Vsub] = RestrictToZ2SymPlusInversion(N, HamC)
%RestrictToZ2SymPlusInversion compute HamC and HamB in symmetric subspace 
%   obeying both Z_2 (X^{\otimes N}) and inversion symmetry
%
%   [HCsym, HBsym, Vsub] = RestrictToZ2SymPlusInversion(N, HamC)
%   
%   HCsym and HBsym are HamC and HamB restricted to symmetric subspace
%       note dimension of symmetric subspace is
%       dimH = 2^(floor(N/2)) + (2^(N-1)-2^(floor(N/2))/2
%
%   Vsub is a 2^N by dimH matrix, listing the orthonormal basis vector of
%   symmetric subspace in the original 2^N-dimension Hilbert space. 
%   e.g. Vsub*Vsub' is the projector onto the symmetric subspace
%

%% work with only first half of the Hilbert space due to X^{\otimes N} symmetry

reduceSym = @(Op) Op(1:end/2,1:end/2) + fliplr(Op(1:end/2,end/2+1:end));


Hsym = reduceSym(HamC);

%% use inversion symmetry to further reduce Hilbert space dimension

sx = sparse([0,1; 1,0]);
sy = sparse([0,-1i; 1i, 0]);
sz = sparse(diag([1,-1]));

InvertSym = speye(2^N);
for ind = 1:floor(N/2)
    InvertSym = InvertSym * 1/2*(Ham2LTerm(sx, sx, ind, N-ind+1, N) + Ham2LTerm(sy, sy, ind, N-ind+1, N) + Ham2LTerm(sz, sz, ind, N-ind+1, N) + speye(2^N));
end

InvertSymSub = reduceSym(InvertSym);

Ipalindrome = find(diag(InvertSymSub) > 0);
[r,c] = find(triu(InvertSymSub,1));

numPalindrome = length(Ipalindrome);
dimHilb = numPalindrome + length(r); % or numPalindrome + (2^(N-1)-numPalindrome)/2;

Vsub = sparse([Ipalindrome; r; c], ...
    [1:numPalindrome, numPalindrome+1:dimHilb, numPalindrome+1:dimHilb], ...
    [ones(numPalindrome, 1); ones(length(r)*2,1)/sqrt(2)], 2^(N-1), dimHilb);


HamB = krondist(sx, N);
HamBSub = reduceSym(HamB);

HBsym = Vsub'*HamBSub*Vsub;
HCsym = Vsub'*Hsym*Vsub;

% generate the basis vectors in the full Hilbert space
Vsub = [Vsub; flipud(Vsub)]/sqrt(2);
% or equivalently, Vsub = [speye(2^(N-1)); kronrec(sx, N-1)]*Vsub/sqrt(2);

end