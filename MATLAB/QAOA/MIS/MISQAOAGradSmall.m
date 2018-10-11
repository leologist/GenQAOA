function [F_ave, F_grad, psiout] = MISQAOAGradSmall(p,HamC,HamB,HamBV,HamBD,param)
%MISQAOAGradSmall computes the average value of Independent Set found using QAOA
%   as well as the gradient of the objective function, suitable for small
%   size Hamiltonians HamC and HamB
%
% [F_ave, F_grad, psiout] = MISQAOAGradSmall(p,HamC,HamB,HamBV,HamBD,param)
%
%   HamBV and HamBD are pre-computed diagonalization of HamB
%   [HamBC, HamBD] = eig(HamB);
%
% F_ave: negative of the objective function (for minimization)
% F_grad: gradient of F_ave with respect to param
% psiout: output wavefunction
%
% first p-1 params are gammas, second p params are betas
% save the psi_p in memory after each cycle to speed up the computation

% EvolHamB = @expv;
% EvolHamB = @expmv;
EvolHamB = @(t,HamBV,HamBD,psi) HamBV*(exp(-1i*t*HamBD).*(HamBV'*psi));

paramC = [0; param(1:p-1)]; paramB = param(p:2*p-1);
F_grad = zeros(numel(paramC)+numel(paramB),1);

psi_in = sparse(1,1,1,numel(HamC),1);
psi_p = zeros(numel(psi_in),2*p+1); % psi_in, psi_1, psi_2, ... psi_p, psi_Mid, psi_p+1, ... psi_2p-1
psi_p(:,1) = psi_in; clear psi_in; 

%------------- the wave function (psi_p) after each cycle
% psi after each cycle, first p steps
psi_p(:,2) = EvolHamB(paramB(1),HamBV,HamBD,psi_p(:,1)); % after first beta
for ind = 2:p
    psi_temp = exp(-1i*paramC(ind)*HamC).*psi_p(:,ind);
    psi_p(:,ind+1) = EvolHamB(paramB(ind),HamBV,HamBD,psi_temp);
end

% psi after each cycle, after HamC, next p steps
psi_p(:,p+2) = HamC.*psi_p(:,p+1); % psi_Mid
for ind = 1:p-1
    psi_temp = EvolHamB(-paramB(p+1-ind),HamBV,HamBD,psi_p(:,ind+p+1));
    psi_p(:,ind+p+2)= exp(1i*paramC(p+1-ind)*HamC).*psi_temp;
end
%------------- the wave function (psi_p) after each cycle

F_ave = real(psi_p(:,p+1)'*psi_p(:,p+2));

for ind = 2:p
    F_grad(ind) = psi_p(:,ind)'*((1i*HamC).*psi_p(:,2*p+3-ind));    
    F_grad(p+ind) = psi_p(:,ind+1)'*((1i*HamB)*psi_p(:,2*p+2-ind));                            
end
F_grad(p+1) = psi_p(:,2)'*((1i*HamB)*psi_p(:,2*p+1));
F_grad(1) = []; % gamma1 is absent

% Make them negative for fmincon (minimization)
F_grad = -2*real(F_grad); % derivative has two parts 
F_ave = -F_ave;

if nargout >= 3
    psiout = psi_p(:, p+1);
end

end

