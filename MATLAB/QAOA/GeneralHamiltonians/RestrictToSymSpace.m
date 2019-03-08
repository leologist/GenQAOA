function [Vsym, Dsym, Hsym] = RestrictToSymSpace(H, flagSym)
%RestrictToSymSpace compute the symmetric subspace eigenvectors, eigenvalues
%   and effective Hamiltonian in the input N-qubit Hamiltonian that has Z_2
%   symmetry with respect to flipping of all spins: [X^{\otimes N}, H]=0
%
%   [Vsym, Dsym, Hsym] = RestrictToSymSpace(H)
%   Vsym, Hsym are 2^{N-1}-dimensional matrices,
%   and Dsym is 2^{N-1} x 1 vector of eigenvalues
%
%   Explanation: Let P=X^{\otimes (N-1)}
%       If state |psi> is Z_2 symmetric, i.e. |psi> = [v; Pv]
%       Then we can restrict to |psi'> = \sqrt{2} v (note it's normalized)
%            and computing e^{-iHt}|psi> is the same as computing
%                     Vsym*e^{-i*Dsym*t}*Vsym'*|psi'>
%            and <a|H|b> = <a'|Hsym|b'>
%
%   More specifically, if H = [H1, H2; H3, H4], then Hsym = H1 + H2*P
%              (Note H1 = P*H4*P, and H2 =P*H3*P if [X^{\otimes N}, H] = 0)



if nargin < 2
    flagSym = true;
end


if flagSym % work in symmetric subspace eig(X^{\otimes N}) = +1
    Hsym = H(1:end/2, 1:end/2) + fliplr(H(1:end/2, end/2+1:end));
else % work in asymmetric subspace eig(X^{\otimes N}) = -1
    Hsym = H(1:end/2, 1:end/2) - fliplr(H(1:end/2, end/2+1:end));
end
[Vsym, Dsym] = eig(full(Hsym));
Dsym = diag(Dsym);


end
