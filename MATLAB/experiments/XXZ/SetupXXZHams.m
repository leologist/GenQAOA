function [QAOAhelperfcn, HamObj, ...
    HamC, HamZ, HamB, EvolC, EvolZ, EvolB] =  SetupXXZHams(N, alphaLR, flagAlternateZ, Delta)
%SetupXXZHams sets up the Hamiltonians and evolution functions for QAOA
%           in one-dimension using ZZ, Z, and X-type driver to optimize the
%           XXZ Hamiltonian
%
%   Arguments: 
%       N = system size
%       alphaLR = power of long-range interaction J_ij = 1/|i-j|^alpha 
%                 (Inf if nearest neighbor only)
%       flagAlternateZ = true if using \sum_i (-1)^i Z_i as driver,
%                        false if using \sum_i Z_i as driver
%       Delta = anisotropy in the XXZ objective Hamiltonian (default=-3/4)
%
%   Returns:
%       QAOAhelperfcn = @(p, param) an evolution function based on
%                       MultiQAOAGrad, returns <HamObj> and its gradient
%                       for a given p and a set of parameters param
%       HamObj = the XXZ objective Hamiltonian, which is
%                \sum_{i=1}^{n-1} (XX + YY + Delta*ZZ)_{i,i+1}
%       HamC = ZZ-type driver (potentially has long-range interaction)
%       HamZ = Z-type driver
%       HamB = X-type driver (\sum_i X_i)
%       EvolC = @(psi, gamma) evolves a given psi by exp(-1i*gamma*HamC)
%       EvolZ = @(psi, alpha) evolves a given psi by exp(-1i*alpha*HamZ)
%       EvolB = @(psi, beta) evolves a given psi by exp(-1i*beta*HamB)
%

if nargin < 4
    Delta = -3/4;
end

%% Set up Hamiltonian
sx = sparse([0,1; 1,0]);
sy = sparse([0,-1i; 1i, 0]);
sz = sparse(diag([1,-1]));

% objective Hamiltonian
HamObj = sparse(0);

for ind = 1:N-1
    HamObj = HamObj + Ham2LTerm(sx, sx, ind, ind+1, N) +...
        Ham2LTerm(sy, sy, ind, ind+1, N) + Delta*Ham2LTerm(sz, sz, ind, ind+1, N);
end

% ZZ driver

HamC = sparse(0);

if alphaLR < Inf
    Jij = @(x, y) 1/abs(x-y)^alphaLR;

    for ind = 1:N-1
        for ind2 = ind+1:N
            HamC = HamC + Jij(ind, ind2) * Ham2LTerm(sz, sz, ind, ind2, N);
        end
    end
else
    for ind = 1:N-1
        HamC = HamC + Ham2LTerm(sz, sz, ind, ind+1, N);
    end
end

if flagAlternateZ
    HamZ = krondist(sz, N, (-1).^(0:N-1));
else
    HamZ = krondist(sz, N);
end


HamB = krondist(sx, N);
HamZdiag = full(diag(HamZ));
HamCdiag = full(diag(HamC));

EvolB = @(psi, beta) EvolHamB(N, beta, psi, false);
EvolC = @(psi, gamma) exp(-1i*gamma*HamCdiag).*psi;
EvolZ = @(psi, alpha) exp(-1i*alpha*HamZdiag).*psi;

psi0 = ones(2^N,1)/sqrt(2^N);

QAOAhelperfcn = @(p, param) MultiQAOAGrad(p, HamObj, {HamC, HamZ, HamB}, param, psi0, {EvolC, EvolZ, EvolB});
