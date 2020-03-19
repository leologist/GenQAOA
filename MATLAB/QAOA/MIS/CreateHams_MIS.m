function [HamC, HamB, stringsLegal] = CreateHams_MIS(N, wGraph, Omega)
%CreateHams_MIS is a helper function to create HamC (objective function)
% and HamB (mixing Hamiltonian) in the constraint subspace of legal strings
% (independent sets) of a graph on N nodes specified by edge list wGraph
%
% Usage:
%  [HamC, HamB] = CreateHams_MIS(N, wGraph)
%  [HamC, HamB] = CreateHams_MIS(N, wGraph, Omega)
%  [HamC, HamB, stringsLegal] = CreateHams_MIS(N, wGraph, Omega)
%
% Output:
%   HamC = a vector, corrseponding to diagonal of the Hamiltonian \sum_i n_i
%
%   HamB = a sparse Hamiltonian matrix from \sum_i Omega_i*X_i
%                    projected in the constraint subspace
%
%   stringsLegal = a table of all independent sets listed by row
%
% Input:
%   N = number of nodes
%   wGraph = list of edges in the graph, give by each row (i, j)
%   Omega = array of individual coupling strengths (default Omega = ones(N,1))
%

if nargin == 2    
    Omega = ones(N,1); % homogeneous onsite strength
end

[stringsLegal, vecLegality] = CreateSubspace(N, wGraph);
HamC = sum(stringsLegal,2); % \sum_i n_i

HamB = CreateHamB(stringsLegal, vecLegality, Omega);

end
