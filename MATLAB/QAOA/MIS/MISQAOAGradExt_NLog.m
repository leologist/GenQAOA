function [G_ave, G_grad, psiout] = MISQAOAGradExt_NLog(p, Obj, HamC, HamB, param)
%MISQAOAGradExt_log computes the -log of the average Obj found using QAOA 
%   evolving under HamC and HamB. It also gives the gradient of -log<HamObj>.
%
%   [G_ave, G_grad, psiout] = MISQAOAGradExt_NLog(p, Obj, HamC, HamB, param)
%
%   Obj is a vector corresponds to objective function in the Z basis
%   HamC is a vector corresponding to the driver in the Z basis
%   HamB is the mixing Hamiltonian
%
% G_ave: -log(<Obj>) (for minimization)
% G_grad: gradient of G_ave with respect to param
% psiout: output wavefunction
%
% first p-1 params are gammas, second p params are betas
% save the psi_p in memory after each cycle to speed up the computation

EvolHamB = @expmv;

param = param(:);

paramC = [0; param(1:p-1)];
paramB = param(p:2*p-1);
G_grad = zeros(numel(paramC)+numel(paramB),1);

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

% psi after each cycle, after HamC, next p steps
psi_p(:,p+2) = Obj.*psi_p(:,p+1); % psi_Mid
for ind = 1:p-1
    psi_temp = EvolHamB(-1i*-paramB(p+1-ind),HamB,psi_p(:,ind+p+1));
    psi_p(:,ind+p+2)= exp(1i*paramC(p+1-ind)*HamC).*psi_temp;
end
%------------- the wave function (psi_p) after each cycle

G_ave = real(psi_p(:,p+1)'*psi_p(:,p+2));

for ind = 2:p
    G_grad(ind) = psi_p(:,ind)'*((1i*HamC).*psi_p(:,2*p+3-ind));    
    G_grad(p+ind) = psi_p(:,ind+1)'*((1i*HamB)*psi_p(:,2*p+2-ind));                            
end
G_grad(p+1) = psi_p(:,2)'*((1i*HamB)*psi_p(:,2*p+1));
G_grad(1) = []; % gamma1 is absent

% Make them negative for fmincon (minimization)
G_grad = -2*real(G_grad)/G_ave; % derivative has two parts 
G_ave = -log(G_ave);

if nargout >= 3
    psiout = psi_p(:, p+1);
end

end

