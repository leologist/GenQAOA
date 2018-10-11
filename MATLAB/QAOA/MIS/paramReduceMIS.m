function paramOut = paramReduceMIS(param,p)
%paramReduceMIS reduces the parameters into the standard range for MIS problem
%
% Symmetries: 1. gamma to -gamma; 2. beta to -beta; 3. gamma_i = gamma_i + pi and beta_j = - beta_j for j>=i  
% So take the gauge: gamma_2 > 0, beta > 0, gamma = -pi to pi
% Input: param = (gamma_2, gamma_3, ..., gamma_p, beta_1, beta_2, ... beta_p); gamma_1 is absent
% p = depth 
% ----- NOTE: this can still has some problem when the parameters are stuck on the boundary for p = 2.

param = param(:);
gamma = param(1:(p-1));
beta = param(p:end);

gamma = mod(gamma,2*pi); 
gamma(gamma>pi) = gamma(gamma>pi) - 2*pi; % range: -pi to pi
gamma = [NaN; gamma(:)]; % add a placeholder for gamma_1
gammaTh = 10^-3;

beta = sign(beta(1))*beta; % beta_1 > 0 

if p>1
    for indp = 2:p
        if beta(indp) < 0
            % shift gamma_i by pi to change the sign for all beta_j for j>=i
            if gamma(indp) < 0
                gamma(indp) = gamma(indp)+pi;
            else
                gamma(indp) = gamma(indp)-pi;
            end            
            beta(indp:p) = -beta(indp:p);             
        end
    end    
end

if p>1
    gamma = sign(gamma(2))*gamma; % gamma_2 > 0
end

if p == 2    
    % if one of the beta is zero, gamma_2 has no use
    if beta(1) < gammaTh || beta(2) < gammaTh
        gamma(2) = 0;
    end    
    % if gamma_2 = 0, or pi, -pi, reduce two betas into one 
    if abs(gamma(2)-pi) < gammaTh || abs(gamma(2)) < gammaTh 
        gamma(2) = 0;
        beta(1) = abs(beta(1) - beta(2));
        beta(2) = 0;            
    end
end

paramOut = [gamma(2:end); beta];

end