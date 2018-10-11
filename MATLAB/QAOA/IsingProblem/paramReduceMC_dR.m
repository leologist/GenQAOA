function paramOut = paramReduceMC_dR(param,d)
% paramReduceMC_dR reduces the QAOA MaxCut parameters into the standard
% range for unweighted d-regular graphs
%   
%   0 < gamma < pi, -pi/4 < beta < pi/4

gamma = param(1:end/2);
beta = param(end/2+1:end);

gamma = mod(gamma,2*pi);
beta = mod(beta, pi/2);

% time reversal symmetry 
if (pi/2 < gamma(1) && gamma(1) < pi) || (3*pi/2 < gamma(1) && gamma(1) < 2*pi)
    gamma = mod(-gamma,2*pi);
    beta = mod(-beta,pi/2);    
end

if mod(d,2) == 0 % even regular
    gamma = mod(gamma,pi);
elseif mod(d,2) == 1 % odd regular 
    [gamma,beta] = paramReduceOdd(gamma,beta);
end


beta = mod(beta, pi/2);
Ilarge = (beta > pi/4) ;
beta(Ilarge) = beta(Ilarge) - pi/2;

paramOut = [gamma(:); beta(:)];

end

% -------------- odd regular graphs -------------- %
function [gamma,beta] = paramReduceOdd(gamma,beta)

for ind = 1:numel(gamma)
    if gamma(ind) > pi
        gamma(ind) = gamma(ind) - pi; 
        beta(ind:end) = pi/2 - beta(ind:end);                
    end
end

end