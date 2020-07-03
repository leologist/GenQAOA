function xy = rand2Dgrid(Ls, fillingFraction, lattice_function)
%rand2Dgrid generates a set of xy coordinates of nodes in a Lx-by-Ly
%   lattice  with a specified filling fraction
%
% Usage:
%   xy = rand2Dgrid(Ls, fillingFraction)
%
% Input: 
%   Ls = either a tuple [Lx, Ly] or a scalar L (then Lx = Ly = L)
%        indicating the linear size of the square lattice
%   fillingFraction = a number between 0 and 1 indicating the percentage of
%        nodes that are filled (default = 1);
%   lattice_function = (optional) a function @(x, y) that takes two integer
%          (x,y) and output the positioon the vertex 
%          (default: lattice_function = @(x,y) [x,y], i.e. square lattice)
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

if nargin < 3
    lattice_function = @(x,y) [x,y];
end

if nargin < 2
    fillingFraction = 1; 
end
    

xy = zeros(Lx*Ly, 2);

for x = 0:Lx-1
    for y = 0:Ly-1
        xy(x*Ly + y + 1, :) = lattice_function(x,y);
    end
end

if fillingFraction < 1
    N = Lx*Ly;
    xy = xy(randperm(N, round(N*fillingFraction)), :);
end



end
