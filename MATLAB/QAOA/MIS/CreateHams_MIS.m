function [HamC, HamB] = CreateHams_MIS(N, wGraph, Omega)
%CreateHams_MIS helper function to create HamC (objective function) and
%   HamB (mixing Hamiltonian) in the constraint subspace of legal strings
%   (independent sets)
%
% [HamC, HamB] = CreateHams_MIS(N, wGraph, Omega)
%   HamC is a vector (corrsponding to diagonal of the Hamiltonian)
%   HamB is a full Hamiltonian matrix from \sum_i Omega_i*s_i^x
%

if nargin == 2    
    Omega = ones(N,1); % homogeneous onsite strength
end

[stringsLegal,vecLegality] = CreateSubspace(N, wGraph);
HamC = sum(stringsLegal,2); % \sum_i n_i

HamB = CreateHamB(stringsLegal, vecLegality, Omega);
end
