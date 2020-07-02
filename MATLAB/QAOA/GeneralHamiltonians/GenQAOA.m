function [F, psi] = GenQAOA(p, HamC, param, psi, EvolC, EvolB)
%GenQAOA performs QAOA circuit with general Hamiltonian
%
%   [F, psi] = GenQAOA(p, HamC, param, psi, EvolC, EvolB)
%       N = number of particles
%       p = QAOA level/depth
%       HamC = problem Hamiltonian (trying to minimize) - a matrix
%       EvolC = @(psi, gamma) ...
%       EvolB = @(psi, beta) ...
%       param = 2*p parameter values, i.e. [gammas, betas]
%       psi = initial state (default = all one vector)
%
%   F = <psi|HamC|psi>
%   psi = output state of QAOA circuit

gammas = param(1:p);
betas = param(p+1:end);

if isempty(psi)
    psi = 1/sqrt(length(HamC))*ones(length(HamC),1);
end
    
for ind = 1:p
    psi = EvolB(EvolC(psi, gammas(ind)), betas(ind));
end

F = real(psi'*HamC*psi);


end