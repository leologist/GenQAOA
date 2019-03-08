function [F, F_grad] = GenQAOAGradSmall(N, p, HamC, HamC_V, HamC_D, param, flagSym)
%GenQAOAGradSmall performs QAOA circuit with general Hamiltonian, and computes
%gradient, for small enough Hamiltonian that can be diagonalized
%
%   [F, F_grad] = GenQAOAGradSmall(N, p, HamC, HamC_V, HamC_D, param, flagSym)
%       N = number of particles
%       p = QAOA level/depth
%       HamC = problem Hamiltonian (trying to minimize)
%       HamC_V, HamC_D are eigenvectors and eigenvalues of HamC
%       param = 2*p parameter values, i.e. [gammas, betas]
%       flagSym = true if using only first half of Hilbert space
%
%   F = <psi|HamC|psi>
%   F_grad = gradient of F with respect to 2*p parameters
%

gammas = param(1:p);
betas = param(p+1:end);
F_grad = zeros(numel(param),1);

if nargin < 7  % Z2 symmetry flag by default is false
    flagSym = false;
end

if flagSym
    psi_in = 1/sqrt(2^(N-1))*ones(2^(N-1),1);
else
    psi_in = 1/sqrt(2^(N))*ones(2^(N),1);
end

EvolC = @(psi, gamma) HamC_V*(exp(-1i*gamma*HamC_D).*(HamC_V'*psi));

EvolB = @(psi, beta) EvolHamB(N, beta, psi, true);

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
    F_grad(ind) = psi_p(:,ind)'*((1i*HamC)*psi_p(:,2*p+3-ind));    
    psi_temp = 0;
    
    if ~flagSym   % Z2 symmetry flag
        for ind2 = 1:N
            psi_temp = psi_temp + 1i*MultSingleSpin(psi_p(:,2*p+2-ind),N,ind2,1); % gradient of exp(1i*HamB*beta)
        end
    else
        for ind2 = 1:(N-1)
            psi_temp = psi_temp + 1i*MultSingleSpin(psi_p(:,2*p+2-ind),N-1,ind2,1); % gradient of exp(1i*HamB*beta)
        end
        psi_temp = psi_temp + 1i*flip(psi_p(:,2*p+2-ind));
    end        
    
    F_grad(p+ind) = psi_p(:,ind+1)'*psi_temp; 
end

F_grad = 2*real(F_grad); % derivative has two parts 

end