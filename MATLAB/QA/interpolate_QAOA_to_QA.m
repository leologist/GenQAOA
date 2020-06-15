function [T, fun, x_interp, y_interp] = interpolate_QAOA_to_QA(param)
%interpolate_QAOA_to_QA takes a parameter vector of length 2p and returns
%   an annealing path made from linear interpolation
%
%   Usage: 
%       [T, fun] = interpolate_QAOA_to_QA(param)
%       [T, fun, x_interp, y_interp] = interpolate_QAOA_to_QA(param)
%
%   Input:
%       param = the 2p parameters of QAOA in the form
%               (gamma_1, ..., gamma_p, beta_1, ..., beta_p)
%
%   Output:
%       T = total time
%       fun = a 1-parameter func that takes s in [0,1] to f(s) in [0,1]
%             such that f(t_i) = gamma_i / (|gamma_i| + |beta_i|)
%             where t_i = \sum_{k=1}^i [ |gamma_k| + |beta_k|) -
%                                        (|gamma_i| + |beta_i|) / 2 ]

p = length(param)/2;
param = param(:); % flattern it to make sure it is an 2p x 1 vector

gammas = param(1:p);
betas = param(p+1:2*p);
alphas = gammas./(gammas+betas);

t_interp = cumsum(gammas+betas);

T = t_interp(end);

x_interp = [0; ([0; t_interp(1:end-1)] + t_interp)/2/T; 1];
y_interp = [0; alphas; 1];

fun = @(s) interp1(x_interp, y_interp, s, 'linear');