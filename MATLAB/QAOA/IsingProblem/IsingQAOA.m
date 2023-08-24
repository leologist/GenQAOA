function [F_ave, psi] = IsingQAOA(N,p,HamC,param,flagSym, psi)
%IsingQAOA computes the average value of an Ising objective function HamC
% running the QAOA algorithm
%
%   [F_ave, psi] = IsingQAOA(N, p, HamC, param)
%   [F_ave, psi] = IsingQAOA(N, p, HamC, param, flagSym)
%   [F_ave, psi] = IsingQAOA(N, p, HamC, param, flagSym, psi)
%
%	HamC is a vector corresponds to diagonal of the Hamiltonian
%
%   flagSym = (default false) denotes a flag for considering the Z2 symmetry
%   (fix first spin to be 0, and reduce Hilbert space size to half)
%
%   psi = initial wavefunction (default: uniform superposition)
%
%   First p params are gammas, second p params are betas.
%   In general, 0<= beta < pi/2. For unweighted graph, 0<= gamma < 2*pi

paramC = param(1:p); paramB = param(p+1:end);

if nargin <= 4
    flagSym = false;
end

if nargin <= 5
    if ~flagSym
        psi = 1/sqrt(2^N)*ones(2^N,1);
    else
        psi = 1/sqrt(2^(N-1))*ones(2^(N-1),1);
    end
end

for ind = 1:p
    psi = exp(-1i*paramC(ind)*HamC).*psi;
    psi = EvolHamB(N,paramB(ind),psi,flagSym);    
end
F_ave = real(psi'*(HamC.*psi)); % make it real explicitly 

end