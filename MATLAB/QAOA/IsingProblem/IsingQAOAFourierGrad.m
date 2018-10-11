function [F_ave, F_grad] = IsingQAOAFourierGrad(N, p, HamC, fParam, flagSym)
%IsingQAOAFourierGrad computes the average value of an Ising objective function HamC
%   running the QAOA algorithm, and efficiently calculates its gradient
%   with respect to "Fourier" parameters fParam
%
% [F_ave, F_grad] = IsingQAOAFourierGrad(N, p, HamC, fParam, flagSym)
%
% HamC is a vector corresponds to diagonal of the Hamiltonian
% fParam:  length 2*q
%          first q parameters are u_k, second q parameters are v_k
%          k-th component correspond to sin(k*pi*i/p)
%          q should be <= p  
% flagSym == 1 means we consider Z2 symmetry of Hilbert space (dim=2^N-1)
%           (default: false)
%
%  Uses IsingQAOAGrad.m


%%
    
if nargin <= 4
    flagSym = false;
end
    
p_trunc = length(fParam)/2;

u = fParam(1:p_trunc); u = u(:)';
v = fParam(p_trunc+1:end); v = v(:)';

[si, k_ind] = meshgrid(1:p, 0:(p_trunc-1));

MA = sin((k_ind + 1/2)*pi.*(si-1/2)/p);
MB = cos((k_ind + 1/2)*pi.*(si-1/2)/p);

gamma = u*MA;
beta = v*MB;

if nargout >= 2
    [F_ave, F_grad_raw] = IsingQAOAGrad(N,p,HamC,[gamma, beta],flagSym);
    grad_gamma = F_grad_raw(1:p);
    grad_beta = F_grad_raw(p+1:end);

    F_grad = [MA*grad_gamma; MB*grad_beta];
else
    F_ave = IsingQAOA(N,p,HamC,[gamma, beta],flagSym);
end

end