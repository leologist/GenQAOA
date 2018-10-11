function ncomp = num_connected_comp(AdjMat)
%num_connected_comp returns number of connected component in the graph
%
%   ncomp = num_connected_comp(AdjMat)


L = diag(sum(AdjMat)) - AdjMat;
eigvals = eig(L);
ncomp = nnz(eigvals < 1e-6);