function xy = rand2Dgrid(Ls, fillingFraction)
%rand2Dgrid generates a set of xy coordinates of nodes in a Lx-by-Ly square
% lattice with a specified filling fraction
%
% Usage:
%   xy = rand2Dgrid(Ls, fillingFraction)
%
% Input: 
%   Ls = either a tuple [Lx, Ly] or a scalar L (then Lx = Ly = L)
%        indicating the linear size of the square lattice
%   fillingFraction = a number between 0 and 1 indicating the percentage of
%        nodes that are filled
% 
% Output:
%   xy = K x 2 matrix of the K nodes in the square lattice, where
%           K = round(Lx*Ly*fillingFraction)

if numel(Ls) == 1
    Lx = Ls; Ly = Ls;
elseif numel(Ls) == 2
    Lx = Ls(1); Ly = Ls(2);
else
    error('Invalid input for Ls')
end

if nargin < 2
    fillingFraction = 1; 
end
    

xy = zeros(Lx*Ly, 2);

for x = 0:Lx-1
    for y = 0:Ly-1
        xy(x*Lx + y + 1, :) = [x,y];
    end
end

N = Lx*Ly;
xy = xy(randperm(N, round(N*fillingFraction)), :);



end

