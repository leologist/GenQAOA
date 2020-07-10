function [F_ave, psi_out] = IsingQAOA(N,p,HamC,param,flagSym)
%IsingQAOA computes the average value of an Ising objective function HamC
% running the QAOA algorithm
%
%   [F_ave, psi_out] = IsingQAOA(N, p, HamC, param)
%   [F_ave, psi_out] = IsingQAOA(N, p, HamC, param, flagSym)
%
%	HamC is a vector corresponds to diagonal of the Hamiltonian
%
%   flagSym = (default false) denotes a flag for considering the Z2 symmetry
%   (fix first spin to be 0, and reduce Hilbert space size to half)
%
%   First p params are gammas, second p params are betas.
%   In general, 0<= beta < pi/2. For unweighted graph, 0<= gamma < 2*pi

paramC = param(1:p); paramB = param(p+1:end);

if nargin <= 4 || ~flagSym
    psi_out = 1/sqrt(2^N)*ones(2^N,1);
    for ind = 1:p
        psi_out = exp(-1i*paramC(ind)*HamC).*psi_out;
        psi_out = EvolHamB(N,paramB(ind),psi_out);
    end
else
    psi_out = 1/sqrt(2^(N-1))*ones(2^(N-1),1);
    for ind = 1:p
        psi_out = exp(-1i*paramC(ind)*HamC).*psi_out;
        psi_out = EvolHamB(N,paramB(ind),psi_out,flagSym);    
    end
end

F_ave = real(psi_out'*(HamC.*psi_out)); % make it real explicitly 

end