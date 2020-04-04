function [edges, AdjMat] = RydbergGraphFromXY(xy, alpha)
%RYDBERGGRAPHFROMXY create graph from inputted xy coordinates, whose edge
%   weights are 1/r_ij^alpha or HeavisideTheta[1-r_ij] if alpha=Inf
%   
%   [edges, AdjMat] = RydbergGraphFromXY(xy)
%   [edges, AdjMat] = RydbergGraphFromXY(xy, alpha)
%
%   optional input: alpha (default = 6)

if nargin < 2
    alpha = 6;
end

NNodes = size(xy,1);
AdjMat = zeros(NNodes);

for ind = 1:NNodes
    for ind2 = ind+1:NNodes
        dist = norm(xy(ind,:) - xy(ind2,:));
        if alpha < Inf
            temp = 1/dist^alpha;
        else
            if dist <= 1
                temp = 1;
            else
                temp = 0;
            end
        end
        AdjMat(ind, ind2) = temp;
    end
end

[row, col, weight] = find(AdjMat);
edges = [row, col, weight];

AdjMat = AdjMat + AdjMat';

end

