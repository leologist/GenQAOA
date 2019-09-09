function [paramC, paramB] = ConvertMaxCutParam(param, direction, real_p, p_trunc)
%ConvertMaxCutParam(param, direction, real_p, p_trunc)
%   Convert between fourier and real space parameters for MaxCut QAOA
%   using DST-4/DCT-4
%
%   [gamma, beta] = ConvertMaxCutParam(param, 1, real_p, p_trunc)
%   [u,v] = ConvertMaxCutParam(param, -1, real_p, p_trunc)
%   
%   direction = 1 if going from [u,v] to [gamma, beta]; -1 from [gamma, beta] 
%   real_p = length of gamma, beta (optional)
%   p_trunc = length of u, v (optional, default = length(param)/2)
%

param = param(:)';

if nargin <= 3 || isempty(p_trunc)
    p_trunc = length(param)/2;
end

if nargin <= 2 || isempty(real_p)
    real_p = p_trunc;
end


[si, k_ind] = meshgrid(1:real_p, 0:(p_trunc-1));
MA = sin((k_ind + 1/2)*pi.*(si-1/2)/real_p);
MB = cos((k_ind + 1/2)*pi.*(si-1/2)/real_p);
 
if direction == 1
    u = param(1:p_trunc);
    v = param(p_trunc+1:end);
    paramC = u*MA;
    paramB = v*MB;
else
    gamma = param(1:real_p)';
    beta = param(real_p+1:end)';
    paramC = (pinv(MA*MA')*MA*gamma)';
    paramB = (pinv(MB*MB')*MB*beta)';
end


if nargout <= 1
    paramC = [paramC, paramB];
end

end

