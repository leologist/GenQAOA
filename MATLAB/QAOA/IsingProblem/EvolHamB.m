function psi_out = EvolHamB(N,beta,psi_in,flagSym)
%   psi_out = EvolHamB(N,beta,psi_in)
%   flagSym == 1 denotes a flag for considering the Z2 symmetry 
%   This function computes exp(-1i*beta* sum_i sx_i)*psi_in, making
%   use of the special structure of the evolution operator
%   Each term in the sum commutes and we can break them up, making
%   use of exp(-1i beta sx) = I cos(beta) - i sin(beta) sx

if nargin <= 3
    flagSym = false;
end

if ~flagSym
    for ind = 1:N
            psi_in = cos(beta)*psi_in - 1i*sin(beta)*MultSingleSpin(psi_in,N,ind,1);
    end
else
    % consider Z2 symmetry and fix the first qubit to be 0 
    for ind = 1:N-1
        psi_in = cos(beta)*psi_in - 1i*sin(beta)*MultSingleSpin(psi_in,N-1,ind,1);
    end
    psi_in = cos(beta)*psi_in - 1i*sin(beta)*flip(psi_in);    
end 

psi_out = psi_in;

end