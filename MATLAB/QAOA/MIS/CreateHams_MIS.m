function [HamC, HamB] = CreateHams_MIS(N, wGraph, Omega)
%CreateHams_MIS is a helper function to create HamC (objective function)
%  and HamB (mixing Hamiltonian) in the constraint subspace of legal strings
%  (independent sets) of a graph on N nodes specified by edge list wGraph
%
%  Usage:
%  [HamC, HamB] = CreateHams_MIS(N, wGraph)
%  [HamC, HamB] = CreateHams_MIS(N, wGraph, Omega)
%
%   HamC is a vector, corrseponding to diagonal of the Hamiltonian \sum_i n_i
%
%   HamB is a sparse Hamiltonian matrix from \sum_i Omega_i*X_i
%                    projected in the constraint subspace
%
%   Omega is an optional input (assumed to all 1s if not supplied)
%

if nargin == 2    
    Omega = ones(N,1); % homogeneous onsite strength
end

[stringsLegal,vecLegality] = CreateSubspace(N, wGraph);
HamC = sum(stringsLegal,2); % \sum_i n_i

HamB = CreateHamB(stringsLegal, vecLegality, Omega);

end
