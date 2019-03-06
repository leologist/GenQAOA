function [F, psi] = GenQAOA(N, p, HamC, HamC_V, HamC_D, param, flagSym)
%GenQAOA performs QAOA circuit with general Hamiltonian
%
%   [F, psi] = GenQAOA(N, p, HamC, HamC_V, HamC_D, param, flagSym)
%       N = number of particles
%       p = QAOA level/depth
%       HamC = problem Hamiltonian (trying to minimize)
%       HamC_V, HamC_D are eigenvectors and eigenvalues of HamC
%       param = 2*p parameter values, i.e. [gammas, betas]
%       flagSym = true if using only first half of Hilbert space
%
%   F = <psi|HamC|psi>
%   psi = output state of QAOA circuit

gammas = param(1:p);
betas = param(p+1:end);

if nargin <= 4  % Z2 symmetry flag by default is false
    flagSym = false;
end

if flagSym
    psi = 1/sqrt(2^(N-1))*ones(2^(N-1),1);
else
    psi = 1/sqrt(2^(N))*ones(2^(N),1);
end
    
for ind = 1:p
    psi = EvolHamB(N, betas(ind), HamC_V*(exp(-1i*gammas(ind)*HamC_D).*(HamC_V'*psi)), flagSym);
end

F = real(psi'*HamC*psi);


end