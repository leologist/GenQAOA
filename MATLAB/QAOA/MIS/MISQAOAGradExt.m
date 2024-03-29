function [F_ave, F_grad, psiout] = MISQAOAGradExt(p, HamObj, HamC, HamB, param)
%MISQAOAGradExt computes the average value of HamObj found using QAOA 
%   evolving under HamC and HamB. It also gives the gradient.
%
%   [F_ave, F_grad, psiout] = MISQAOAGradExt(p,HamObj,HamC,HamB,param)
%   HamObj is the objective Hamiltonian (for minimizing) which can be a
%           vector or a matrix
%   HamC is a vector corresponding to the driver in the Z basis
%   HamB is the mixing Hamiltonian in the subspace of independent set
%       states generated from sigma_x
%
% F_ave: the objective function
% F_grad: gradient of F_ave with respect to param
% psiout: output wavefunction
%
% first p-1 params are gammas, second p params are betas
% save the psi_p in memory after each cycle to speed up the computation

EvolHamB = @expmv;

paramC = [0; param(1:p-1)];
paramB = param(p:2*p-1);
F_grad = zeros(numel(paramC)+numel(paramB),1);

psi_in = sparse(1,1,1,numel(HamC),1);
psi_p = zeros(numel(psi_in),2*p+1); % psi_in, psi_1, psi_2, ... psi_p, psi_Mid, psi_p+1, ... psi_2p-1
psi_p(:,1) = psi_in; clear psi_in; 

%------------- the wave function (psi_p) after each cycle
% psi after each cycle, first p steps
psi_p(:,2) = EvolHamB(-1i*paramB(1),HamB,psi_p(:,1)); % after first beta
for ind = 2:p
    psi_temp = exp(-1i*paramC(ind)*HamC).*psi_p(:,ind);
    psi_p(:,ind+1) = EvolHamB(-1i*paramB(ind),HamB,psi_temp);
end

% Apply HamObj
if size(HamObj, 2) == 1
    psi_p(:,p+2) = HamObj.*psi_p(:,p+1); % psi_Mid
else
    psi_p(:,p+2) = HamObj*psi_p(:,p+1); % psi_Mid
end

% psi after each cycle, after HamObj, next p steps
for ind = 1:p-1
    psi_temp = EvolHamB(-1i*-paramB(p+1-ind),HamB,psi_p(:,ind+p+1));
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

F_grad = 2*real(F_grad); % derivative has two parts 

if nargout >= 3
    psiout = psi_p(:, p+1);
end

end

