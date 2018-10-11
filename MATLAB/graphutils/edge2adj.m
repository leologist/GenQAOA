function AdjMat = edge2adj(edges, NNodes)
%EDGE2ADJ convert edge list to adjacency matrix
%   AdjMat = edge2adj(edges)

if ~exist('NNodes', 'var')
    NNodes = max(max(edges(:,1)), max(edges(:,2)));
end

if size(edges, 2) < 3
    edges = [edges, ones(size(edges,1),1)];
end

AdjMat = sparse(edges(:,1), edges(:,2), edges(:, 3), NNodes, NNodes);
AdjMat = AdjMat + AdjMat.';

end
