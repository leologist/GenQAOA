function [Xcut, Ycut, XNcut, YNcut, X_A, Y_A, X_B, Y_B] = gplotcut(A, y, coords)
%GPLOTCUT plot graph and a cut specified by vertex labelling
%   GPLOTCUT(A, y, coords)
%   [Xcut, Ycut, XNcut, YNcut, X_A, Y_A, X_B, Y_B] = GPLOTCUT(A, y, coords)
%
%   Input:
%       A = adjacency matrix, n x n
%       y = vertex labeling, n x 1, (+1 = group A, -1 = group B)
%       coords = [x, y] coordinates of vertices, n x 2
%                (optional, default = unit circle of points)
%
%   Output:
%       Xcut, Ycut = NaN punctuated coordinates of edges cut
%       XNcut, YNcut = NaN punctuated coordinates of edges not cut
%       X_A, Y_A = coordinates of vertices in group A
%       X_B, Y_B = coordinates of vertices in group B
%
%   author: Leo Zhou
%   date: 06/28/2017

n = length(A);

if nargin < 3 || isempty(coords)
    coords = [cos((0:n-1)*2*pi/n); sin((0:n-1)*2*pi/n)]';
end

if nargin < 2 || isempty(y)
    y = ones(n,1);
end

[row, col] = find(triu(A));
boolCut = (y(row) .* y(col) ~= 1);

rCut = row(boolCut);
cCut = col(boolCut);
rNCut = row(~boolCut);
cNCut = col(~boolCut);

Xcut = [coords(rCut,1), coords(cCut, 1)];
Ycut = [coords(rCut,2), coords(cCut, 2)];
XNcut = [coords(rNCut,1), coords(cNCut, 1)];
YNcut = [coords(rNCut,2), coords(cNCut, 2)];

Xcut = [Xcut'; NaN(1, length(rCut))];
Ycut = [Ycut'; NaN(1, length(rCut))];
XNcut = [XNcut'; NaN(1, length(rNCut))];
YNcut = [YNcut'; NaN(1, length(rNCut))];

Xcut = Xcut(:); Ycut = Ycut(:);
XNcut = XNcut(:); YNcut = YNcut(:);

ind_A = find(y==1);
ind_B = find(y~=1);
X_A = coords(ind_A, 1);
Y_A = coords(ind_A, 2);
X_B = coords(ind_B, 1);
Y_B = coords(ind_B, 2);

plot(Xcut, Ycut, '-r', XNcut, YNcut, '-b', X_A, Y_A, 'ob')
hold on
plot(X_B, Y_B, '.r','markersize',16);
hold off


end

