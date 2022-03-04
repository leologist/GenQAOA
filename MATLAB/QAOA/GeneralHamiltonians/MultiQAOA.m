function [F, psi] = MultiQAOA(p, HamObj, params, psi, Evolfuns)
%MultiQAOA performs QAOA with multiple evolution generators
%
%   [F, psi] = MultiQAOA(p, HamObj, params, psi_in, Evolfuns)
%       p = QAOA level/depth
%       HamObj = objective Hamiltonian (trying to minimize)
%       EvolHams = cell array of K Hamiltonians used as evolution generator
%       param = p x K matrix parameter values
%       psi_in = initial wavefunction (default or [] = all one vector)
%       Evolfuns = {@(psi, angle) ..., ...} (default: use Krylov expmv)
%                   cell array of function handles
%
%   F = <psi|HamObj|psi>
%   psi = output wavefunction of the circuit

if isempty(psi)
    psi = 1/sqrt(length(HamObj))*ones(length(HamObj),1);
end


numGen = length(Evolfuns);

if numel(params) ~= p*numGen
    error('MultiQAOA: Invalid # params')
end

for q = 1:p
    for k = 1:numGen
        psi = Evolfuns{k}(psi, params(q, k));
    end
end

F = real(psi'*HamObj*psi);



end