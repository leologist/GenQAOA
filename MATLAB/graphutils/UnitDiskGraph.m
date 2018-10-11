function [edges, AdjMat, xy] = UnitDiskGraph(varargin)
%UnitDiskGraph creates unit disk graph
%
%   Usages:
%   [edges, AdjMat, xy] = UnitDiskGraph(xy)
%       generates unit disk graph from coordinates of vertices, where any
%       two vertices with <= 1 distance are connected by an edge
%
%   [edges, AdjMat, xy] = UnitDiskGraph(xy, rho, 1)
%       generates unit disk graph from coordinates of vertices, where any
%       two vertices with <= 1 distance on a LxL torus [L=sqrt(N/rho)] are
%       connected by an edge [N = size(xy,1)]
%
%   [edges, AdjMat, xy] = UnitDiskGraph(N, rho)
%       generates unit disk graph of N vertices with target density rho (open boundaries) 
% 
%   [edges, AdjMat, xy] = UnitDiskGraph(N, rho,0)
%       generates unit disk graph of N vertices with target density rho (open boundaries)
%
%   [edges, AdjMat, xy] = UnitDiskGraph(N, rho,1)
%       generates unit disk graph of N vertices with target density rho (periodic boundaries)
%
%   Output: edges = list of edges (i, j)
%           AdjMat = adjacency matrix
%           xy = coordinates of vertices


if numel(varargin{1}) == 1
    NNodes = varargin{1};
    rho = varargin{2};
    L = sqrt(NNodes/rho);
    xy = rand(NNodes, 2)*L;
else
    xy = varargin{1};
    NNodes = size(xy,1);
    
    if nargin >= 2
        rho = varargin{2};
        L = sqrt(NNodes/rho);
    end
end

AdjMat = zeros(NNodes);

if nargin <= 2 || (nargin == 3 && varargin{3} == 0) % open boundary condition
    for ind = 1:NNodes
        myXY = repmat(xy(ind, :), NNodes, 1);
        dist = sum((xy - myXY).^2,2);
        AdjMat(dist <= 1, ind) = 1;
        AdjMat(ind, ind) = 0;
    end
end

if nargin == 3 && varargin{3} == 1 % periodic boundary condition
    for ind = 1:NNodes
        myXY = repmat(xy(ind, :), NNodes, 1);
        dist = mod(xy-myXY,L);
        dist = min(dist,L-dist); % minimum distance
        dist = sum(dist.^2,2);
        AdjMat(dist <= 1, ind) = 1;
        AdjMat(ind, ind) = 0;
    end
end

[r, c] = find(triu(AdjMat));
edges = [r, c];