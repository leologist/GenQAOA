function paramOut = paramReduceWMC(param)
%paramReduceWMC reduces the QAOA parameters for MaxCut on weighted graphs
%
%   Usage:  paramOut = paramReduceWMC(param)
%   
%   Input format: param = [gamma; beta];
%
%   Output format: paramOut = [gammaOut; betaOut];
%       where gammaOut(end) > 0 and -pi/4 < betaOut < pi/4

param = param(:);
p = length(param)/2;

gamma = param(1:p);
beta = param(p+1:end);

if gamma(end) < 0
    gamma = -gamma;
    beta = -beta;
end

beta = mod(beta, pi/2);
Ilarge = (beta > pi/4) ;
beta(Ilarge) = beta(Ilarge) - pi/2;

paramOut = [gamma; beta];

end

