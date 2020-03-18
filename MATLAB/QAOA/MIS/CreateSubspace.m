function [stringsLegal,vecLegality] = CreateSubspace(N, wGraph)
%CreateSubspace  generates the subspace of legal strings (corresponding to
% independent sets); when two vertex has an edge, they can't be both in 1
% (Rydberg blockade).
%
% It uses the Bron-Kerbosch algorithm to find all maximal independent set.
% It treats all edges as same weight.
%
%   [stringsLegal,vecLegality] = CreateSubspace(N, wGraph)
%
% Input:
%       N = the number of vertex.
%       wGraph = edge lists, in the format of (i,j) or (i,j,w_ij) as rows.
% Output: 
%       stringsLegal = a KxN table indicating all independent sets (legal
%                     strings)
%
%       vecLegality = 2^N x 1 vector, where nonzero numbers occupy the 
%            independent set states in the original Hilbert space, each 
%            labelling the row position of the legal string in the order 
%            specified in the stringsLegal table.
%
% WARNING: THIS FUNCTION MAY NOT WORK AS INTENDED WHEN N>=53, since
%          double has only ~53 bit resolution, so 2^53+1 = 2^53
%          Use GetIndepedentSets.m instead

if isempty(wGraph)
    adjMat = zeros(N); % no edges
else
    adjMat = full(sparse(wGraph(:,1),wGraph(:,2),ones(size(wGraph,1),1),N,N)); % works only when wGraph(:,3) = 1
end

adjMat = adjMat + adjMat';

% the adjacency matrix for the complement graph
adjMatComp = 1- adjMat;
adjMatComp = adjMatComp - eye(N);

MaxIS = maximalCliques(adjMatComp); % max clique in the complement graph is the max ind set for the graph

%-------------------- Find all independent sets from Maximum independent sets --------------------%
MaxIS_Size = max(cellfun(@(x) numel(x), MaxIS));
% pre-generating the binary strings for a MaxIS size
binaryString = cell(MaxIS_Size,1);
for indSize = 1:MaxIS_Size
    binaryString{indSize} = d2b(0:2^indSize-1,indSize);
end

stringNoAll = zeros(2^24,1); % the position of the string in the original 2^N order
count = 1;
for indMaxIS = 1:numel(MaxIS)
    MaxOneSet = MaxIS{indMaxIS}; % one maximal independent set
    MaxOneSetSize = numel(MaxOneSet);
    binaryFactor = 2.^(MaxOneSet-1)';
    stringNoT = binaryString{MaxOneSetSize}*binaryFactor + 1; % the position of the string in the original 2^N order
    stringNoTSize = numel(stringNoT);
    stringNoAll(count:count+stringNoTSize-1) = stringNoT;
    count = count + stringNoTSize;
end

stringNoAll(count:end) = [];
stringNoAll = unique(stringNoAll);

stringsLegal = logical(d2b(stringNoAll-1,N));
vecLegality = sparse(stringNoAll,1,1:numel(stringNoAll),2^N,1);

end


function s=d2b(d,n)
% Adapted from DEC2BIN function

d = d(:); % Make sure d is a column vector.

% Actual algorithm
s = rem(floor(d*pow2(0:-1:1-n)),2);
end










