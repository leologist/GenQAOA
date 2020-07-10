function paramOut = paramINTERP_QAOA_MIS(param, p_out)
%%paramINTERP_QAOA_MIS use the interpolation strategy to convert a set of
%   optimal QAOA parameters for MIS to a good guess at a different p
%
%   Usage:
%       paramOut = paramINTERP_QAOA_MIS(param, p_out)
%
%   param = input parameters (contains 2*p-1 numbers)
%   p_out = the desired level p for the output parameters
%   paramOut = a guess for parameters at p=p_out (a vector of length 2*p_out-1)
%
    p = (numel(param)+1)/2;
    
    gamma0 = param(1:p-1);
    beta0 = param(p:end);
    gamma = interp1(linspace(0,1,p-1), gamma0, linspace(0,1,p_out - 1));
    beta = interp1(linspace(0,1,p), beta0, linspace(0,1,p_out));
    paramOut = [gamma(:); beta(:)];
end