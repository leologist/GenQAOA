function [F, F_grad] = GenQAOAGrad(p, HamC, HamB, param, psi_in, EvolC, EvolB)
%GenQAOAGrad performs QAOA circuit with general Hamiltonian, and computes
%gradient
%
%   [F, F_grad] = GenQAOAGrad(p, HamC, HamB, param, psi_in, EvolC, EvolB)
%       N = number of particles
%       p = QAOA level/depth
%       HamC = problem Hamiltonian (trying to minimize)
%       HamB = problem Hamiltonian (trying to minimize)
%       param = 2*p parameter values, i.e. [gammas, betas]
%       psi_in = initial wavefunction (default=all one vector)
%
%       Optional: EvolC and EvolB are anonymous functions for evolving
%       wavefunctions according to HamC and HamB
%       EvolC = @(psi, gamma) ... (default: use Krylov expmv)
%       EvolB = @(psi, beta) ... (default: use Krylov expmv)
%
%   F = <psi|HamC|psi>
%   F_grad = gradient of F with respect to 2*p parameters
%

gammas = param(1:p);
betas = param(p+1:end);
F_grad = zeros(numel(param),1);


if nargin < 5
    psi_in = 1/sqrt(length(HamC))*ones(length(HamC),1);
end

if nargin < 6
    EvolC = @(psi, gamma) expmv(-1i*gamma, HamC, psi);
end

if nargin < 7
    EvolB = @(psi, beta) expmv(-1i*beta, HamB, psi);
end

psi_p = zeros(numel(psi_in), 2*p+2);
psi_p(:,1) = psi_in; clear psi_in;

for ind = 1:p
    psi_p(:,ind+1) = EvolB(EvolC(psi_p(:,ind), gammas(ind)), betas(ind));
end

psi_p(:,p+2) = HamC*psi_p(:, p+1);

for ind = 1:p
    psi_p(:, ind+p+2) = EvolC(EvolB(psi_p(:,ind+p+1), -betas(p+1-ind)), -gammas(p+1-ind));
end

F = real(psi_p(:,p+1)'*psi_p(:,p+2));


for ind = 1:p
    F_grad(ind) = psi_p(:,ind)'*(1i*HamC*psi_p(:,2*p+3-ind));
    F_grad(p+ind) = psi_p(:,ind+1)'*1i*HamB*psi_p(:,2*p+2-ind);
end

F_grad = 2*real(F_grad); % derivative has two parts 

end