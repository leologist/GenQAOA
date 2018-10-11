function [F_ave, psiout] = MISQAOA(p,HamC,HamB,param)
%MISQAOA computes the average value of Independent Set found using QAOA
%
%   [F_ave, psiout] = MISQAOA(p,HamC,HamB,param)
%   HamC is a vector corresponds to objective function in the computational
%       basis
%   HamB is the mixing Hamiltonian in the subspace of independent set
%       states generated from sigma_x
%
% F_ave: negative of the objective function (for minimization)
% psiout: output wavefunction
% 
% first p-1 params are gammas, second p params are betas
%
% This is more memory-efficient to get [F_ave, psiout] than MISQAOAGrad.m

EvolHamB_local = @expmv;

param = param(:);
paramC = [0; param(1:p-1)];
paramB = param(p:2*p-1);

psiout = zeros(numel(HamC),1);
psiout(1) = 1;

psiout = EvolHamB_local(-1i*paramB(1),HamB, psiout);

for ind = 2:p
    psiout = exp(-1i*paramC(ind)*HamC).*psiout;
    psiout = EvolHamB_local(-1i*paramB(ind),HamB,psiout);
end

F_ave = -real(psiout'*(HamC.*psiout));


end

