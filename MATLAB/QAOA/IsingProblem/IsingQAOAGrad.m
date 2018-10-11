function [F_ave, F_grad] = IsingQAOAGrad(N,p,HamC,param,flagSym)
%IsingQAOAGrad computes the average value of an Ising objective function HamC
% running the QAOA algorithm, as well as the gradient with respect to
% parameters
%
% [F_ave, F_grad] = IsingQAOAGrad(N,p,HamC,param,flagSym)
%
% Input:
%   HamC is a vector corresponds to diagonal of the Hamiltonian
%   first p params are gammas, second p params are betas, 0<gamma<2*pi, 0<beta<pi/2
%   flagSym == 1 denotes a flag for considering the Z2 symmetry 
%
% Output:
%   F_ave: negative of the objective function (for fmincon)
%   F_grad: negative of the gradient of F_ave
% save the psi_p in memory after each cycle to speed up the computation
% memory usage ~ 2^N*(2*p+4)*16; N=22, p = 30, ~ 4.3G 

paramC = param(1:p); paramB = param(p+1:2*p);
F_grad = zeros(numel(param),1);

if nargin <= 4  % Z2 symmetry flag
    flagSym = false;
end

if ~flagSym
    psi_in = 1/sqrt(2^N)*ones(2^N,1);
    EvolHamBLocal = @(N,beta,psi_in) EvolHamB(N,beta,psi_in);
else
    psi_in = 1/sqrt(2^(N-1))*ones(2^(N-1),1);
    EvolHamBLocal = @(N,beta,psi_in) EvolHamB(N,beta,psi_in,flagSym);
end

psi_p = zeros(numel(psi_in),2*p+2); % psi_in, psi_1, psi_2, ... psi_p, psi_Mid, psi_p+1, ... psi_2p
psi_p(:,1) = psi_in; clear psi_in; 

%------------- the wave function (psi_p) after each cycle
% psi after each cycle, first p steps
for ind = 1:p
    psi_temp = exp(-1i*paramC(ind)*HamC).*psi_p(:,ind);
    psi_p(:,ind+1) = EvolHamBLocal(N,paramB(ind),psi_temp);
end

% psi after each cycle, after HamC, next p steps
psi_p(:,p+2) = HamC.*psi_p(:,p+1); % psi_Mid
for ind = 1:p
    psi_temp = EvolHamBLocal(N,-paramB(p+1-ind),psi_p(:,ind+p+1));
    psi_p(:,ind+p+2)= exp(1i*paramC(p+1-ind)*HamC).*psi_temp;
end
%------------- the wave function (psi_p) after each cycle

F_ave = real(psi_p(:,p+1)'*psi_p(:,p+2));

for ind = 1:p
    F_grad(ind) = psi_p(:,ind)'*((1i*HamC).*psi_p(:,2*p+3-ind));    
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

% Make them negative for fmincon (minimization)
F_grad = -2*real(F_grad); % derivative has two parts 
F_ave = -F_ave;

end

