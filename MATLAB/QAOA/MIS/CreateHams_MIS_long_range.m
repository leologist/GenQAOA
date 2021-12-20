function [HamC, HamB, stringsLegal] = CreateHams_MIS_long_range(xy, alpha, Omega)
%CreateHams_MIS_long_range is a helper function to create HamC (objective function)
%  and HamB (mixing Hamiltonian) in the constraint subspace of legal strings
%  (independent sets) on a Unit Disk graph
%
%  Usage:
%  [HamC, HamB] = CreateHams_MIS(xy)
%  [HamC, HamB] = CreateHams_MIS(xy, alpha)
%  [HamC, HamB] = CreateHams_MIS(xy, alpha, Omega)
%  [HamC, HamB, stringsLegal] = CreateHams_MIS_long_range(xy, alpha, Omega)
%
%  Output:
%   HamC = a vector, corrseponding to diagonal of the Hamiltonian
%                     \sum_i n_i
%   
%   HamB = a sparse Hamiltonian matrix from \sum_i Omega_i*X_i + \sum_{ij} V_ij n_i n_j
%                    projected in the constraint subspace
%                     where V_ij = 1/(r_ij)^alpha
%
%   stringsLegal = a table of all independent sets listed by row
%
%  Input:
%   xy = positions of nodes for the unit disk graph
%   alpha = power of the decaying interaction (default alpha=6)
%   Omega = array of individual coupling strengths (default Omega = ones(N,1))
%

N = size(xy,1);

if nargin < 3
    Omega = ones(N,1); % homogeneous onsite strength
end

if nargin < 2
    alpha = 6;
end


wGraph = UnitDiskGraph(xy);

stringsLegal = GetIndependentSets(N, wGraph);

%% build objective function Hamiltonian, constraint to IS subspace
HamC = sum(stringsLegal,2); % \sum_i n_i

%% build mixing Hamiltonian

if N < 48
    inds = stringsLegal*2.^(0:N-1)' + 1;
    vecLegality = sparse(inds, 1, 1:size(stringsLegal, 1), 2^N, 1);    
    HamB = CreateHamB(stringsLegal, vecLegality, Omega);
else
    HamB = CreateHamB_largeN(stringsLegal, Omega);
end

%% Adding long-range tail
% pre-computing V_ij for all pairs (i,j) for faster performance
Vij = zeros(N,N);
for ind = 1:N-1
    for ind2 = ind+1:N
        Vij(ind, ind2) = 1/norm(xy(ind,:) - xy(ind2,:))^alpha;
    end
end

% Vij(Vij > 1) = 1; % ignore interactions within unit radius as those states are not accessible
    

for ind = 1:size(stringsLegal, 1)
    x = stringsLegal(ind,:);
%     HamC(ind) = HamC(ind) - x*Vij*x';
    HamB(ind, ind) = HamB(ind, ind) + x*Vij*x';
end
