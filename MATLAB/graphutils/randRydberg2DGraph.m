function [edges, AdjMat, xy] = randRydberg2DGraph(NNodes, minDist)
%RANDRYDBERG2DGRAPH generate random graph based on uniformly distributed
%   Rydberg atoms in 2D box with minimum separation and 1/r^6 interaction
%   
%   [AdjMat, xy] = randRydberg2DGraph(NNodes, minDist)
%   
%   NNodes = number of vertices
%   minDist = minimum separation distance between points
%             (optional, default = 0)

boxLength = sqrt(NNodes);
xy = zeros(NNodes,2);
xy(1,:) = rand(1,2)*boxLength;
NFound = 1;

if nargin < 2
    minDist = 0;
end


minDistSq = minDist^2;

ntrials = 0;
while NFound < NNodes
    thisX = rand*boxLength;
    thisY = rand*boxLength;
    distancesSq = (thisX-xy(1:NFound,1)).^2 + (thisY-xy(1:NFound,2)).^2;
    if min(distancesSq) >= minDistSq
        NFound = NFound + 1;
        xy(NFound, :) = [thisX, thisY];
    end
    ntrials = ntrials + 1;
    if mod(ntrials, 1e5) == 0
        fprintf('randRydberg2DGraph: ntrials = %i\n', ntrials);
        NFound = NFound - 1;
        if ntrials >= 1e7
            error('randRydberg2DGraph: exceeded 1e7 tries, quitting')
        end
    end
end

[edges, AdjMat] = RydbergGraphFromXY(xy);


end

