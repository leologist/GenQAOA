function [F, F_grad] = GenQAOAGrad(N, p, HamC, HamC_V, HamC_D, param, flagSym)
%GenQAOAGrad performs QAOA circuit with general Hamiltonian, and computes
%gradient
%
%   [F, F_grad] = GenQAOAGrad(N, p, HamC, HamC_V, HamC_D, param, flagSym)
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

if nargin <= 6  % Z2 symmetry flag by default is false
    flagSym = false;
end

if flagSym
    psi_in = 1/sqrt(2^(N-1))*ones(2^(N-1),1);
else
    psi_in = 1/sqrt(2^(N))*ones(2^(N),1);
end

EvolHamC = @(psi, gamma) HamC_V*(exp(-1i*gamma*HamC_D).*(HamC_V'*psi));


psi_p = zeros(numel(psi_in), 2*p+2);
psi_p(:,1) = psi_in; clear psi;

for ind = 1:p
    psi_p(:,ind+1) = EvolHamB(N, betas(ind), EvolHamC(psi_p(:,ind), gammas(ind)), flagSym);
end

psi_p(:,p+2) = HamC*psi_p(:, p+1);

for ind = 1:p
    psi_p(:, ind+p+2) = EvolHamC(EvolHamB(N, -betas(p+1-ind), psi_p(:,ind+p+1), flagSym), -gammas(p+1-ind));
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