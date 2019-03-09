function [F, F_grad, psi] = GenQAOAGrad(p, HamC, HamB, param, psi_in, EvolC, EvolB)
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

psi = zeros(numel(psi_in), 2*p+2);
psi(:,1) = psi_in; clear psi_in;

for ind = 1:p
    psi(:,ind+1) = EvolB(EvolC(psi(:,ind), gammas(ind)), betas(ind));
end

psi(:,p+2) = HamC*psi(:, p+1);

for ind = 1:p
    psi(:, ind+p+2) = EvolC(EvolB(psi(:,ind+p+1), -betas(p+1-ind)), -gammas(p+1-ind));
end

F = real(psi(:,p+1)'*psi(:,p+2));


for ind = 1:p
    F_grad(ind) = psi(:,ind)'*(1i*HamC*psi(:,2*p+3-ind));
    F_grad(p+ind) = psi(:,ind+1)'*1i*HamB*psi(:,2*p+2-ind);
end

F_grad = 2*real(F_grad); % derivative has two parts

if nargout >= 3
    psi = psi(:,p+1); % output the state at the end of evolution
end

end