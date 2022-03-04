function [F, F_grad, psi] = MultiQAOAGrad(p, HamObj, EvolHams, params, psi_in, Evolfuns)
%MultiQAOAGrad performs multi-QAOA circuit with general objective function,
%   and arbitrary number of evolution generators
%
%   [F, F_grad, psi] = MultiQAOAGrad(p, HamObj, EvolHams, params, psi_in, Evolfuns)
%       p = QAOA level/depth
%       HamObj = objective Hamiltonian (trying to minimize)
%       EvolHams = cell array of K Hamiltonians used as evolution generator
%       param = p x K matrix parameter values
%       psi_in = initial wavefunction (default or [] = all one vector)
%       Evolfuns = {@(psi, angle) ..., ...} (default: use Krylov expmv)
%                   cell array of function handles, same size as EvolHams 
%               used if there's an implementation more efficient than expmv
%
%   F = <psi|HamObj|psi>
%   F_grad = gradient of F with respect to p*K parameters (p x K matrix)
%   psi = output wavefunction of the circuit
%
%   memory usage: (2*(K-1)*p+1) * numel(psi_in) doubles for storaging
%                 wavefunction copies

numGen = numel(EvolHams); % K

if numel(params) ~= p*numGen
    error('Invalid input - not enough parameters for the specified p=%d', p);
end


if nargin < 5
    psi_in = 1/sqrt(length(HamC))*ones(length(HamC),1);
end

if nargin < 6
    Evolfuns = cell(numGen,1);
    for k = 1:numGen
        Evolfuns{k} = @(psi, param) expmv(-1i*param, EvolHams{k}, psi);
    end
end


%% evolving forward
%    uses format: Evolfun{k}(psi, params(q, k))

psiforward = zeros(numel(psi_in), p, numGen-1);

% markers_f = cell(p, numGen-1); % for debugging purpose - commented out for now

for q = 1:p
    
    if q == 1
        psiforward(:, q, 1) = Evolfuns{2}(Evolfuns{1}(psi_in, ...
            params(q, 1)), params(q, 2));
%         markers_f{q, 1} = sprintf('(%d,%d)(%d,%d)', q, 2, q, 1);
    else
        psiforward(:, q, 1) = Evolfuns{2}(Evolfuns{1}(psiforward(:, q-1, end),...
            params(q, 1)), params(q, 2));
%         markers_f{q, 1} = [sprintf('(%d,%d)(%d,%d)', q, 2, q, 1), markers_f{q-1,end}];
    end
    
    for k = 2:numGen-1
        psiforward(:, q, k) = Evolfuns{k+1}(psiforward(:, q, k-1), params(q, k+1));
%         markers_f{q, k} = [sprintf('(%d,%d)', q, k+1), markers_f{q, k-1}];
    end
end



%% evoling backwards after applying HamObj

psibackward = zeros(numel(psi_in), p, numGen-1);

% markers_b = cell(p, numGen-1);

psibackward(:, p, 1) = HamObj*psiforward(:, end, end);
% markers_b{p, 1} = ['[H_O]', markers_f{end,end}];

for k = 2:numGen-1
    psibackward(:, p, k) = Evolfuns{numGen+2-k}(psibackward(:, p, k-1), -params(p, numGen+2-k));
%     markers_b{p, k} = [sprintf('(%d,%d)', p, numGen+2-k), markers_b{p, k-1}];
end

for q = p-1:-1:1
    
    psibackward(:, q, 1) = Evolfuns{1}(Evolfuns{2}(psibackward(:, q+1, end), ...
        -params(q+1, 2)), -params(q+1, 1));
    
%     markers_b{q, 1} = [sprintf('(%d,%d)(%d,%d)', q+1, 2, q+1, 1), markers_b{q+1, end}];
    
    for k = 2:numGen-1
        psibackward(:, q, k) = Evolfuns{numGen+2-k}(psibackward(:, q, k-1), -params(q, numGen+2-k));
%         markers_b{q, k} = [sprintf('(%d,%d)', q, numGen+2-k), markers_b{q, k-1}];
    end
end


%%

F = real(psibackward(:, p, 1)'*psiforward(:, p, numGen-1));

F_grad = zeros(p, numGen);

for q = 1:p
    if q == 1
        F_grad(q, 1) = Evolfuns{1}(Evolfuns{2}(psibackward(:, 1, end), ...
            -params(1, 2)), -params(1, 1))' ...
            * 1i*EvolHams{1}*psi_in;
    else
        F_grad(q, 1) = psibackward(:, q-1, 1)'*1i*EvolHams{1}*psiforward(:, q-1, end);
%         fprintf(2, '%s\n%s\n --- %0.5f\n', markers_f{q-1,end}, markers_b{q-1,1}, F_grad(q, 1))
    end
    
    for k = 2:numGen
        F_grad(q, k) = psibackward(:, q, numGen+1-k)'*1i*EvolHams{k}*psiforward(:, q, k-1);
    end
end

F_grad = -2*real(F_grad); % derivative has two parts

if nargout >= 3
    psi = psiforward(:, end, end); % output the state at the end of evolution
end

end