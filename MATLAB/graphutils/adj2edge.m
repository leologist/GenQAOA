function edges = adj2edge(AdjMat)
%ADJ2EDGE convert adjacency matrix to edges list
%   edges = adj2edge(AdjMat)
%   edges(:,1) is connected to edges(:,2) with weight edges(:,3)

[row, col, v] = find(triu(AdjMat));

edges = [row, col, v];

end
