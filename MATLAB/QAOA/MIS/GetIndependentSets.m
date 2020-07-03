function stringsLegal = GetIndependentSets(N, wGraph)
%GetIndependentSets get all the independent sets from maximal independent
%   sets as a logical array
%
%   Usage:
%   stringsLegal = GetIndependentSets(N, wGraph)
%   
%   Input:
%       N = number of vertices
%       wGraph = list of edges of the graph,
%                whose rows are in the form [i, j] or [i, j, w_ij] 
%
%   Output:            
%      stringsLegal is KxN logical array, K is the number of independent sets
%
%   @ Leo Zhou -- Nov 6, 2018
%

%% Getting all maximal independent set

adjMat = edge2adj(wGraph);

% the adjacency matrix for the complement graph
adjMatComp = ~adjMat - eye(N);

maximalIS = maximalCliques(adjMatComp);

ISsizes = cellfun(@numel, maximalIS);

MISsize = max(ISsizes);


%% pre-generating the binary strings for a MaxIS size
binaryString = cell(MISsize,1);
for IS_size = 1:MISsize
    binaryString{IS_size} = d2b(1:2^IS_size-1,IS_size);
end

%% take each maximal independent set and add all its subset to stringsDecimal


numBits = 50; % needs to be < 53, since MATLAB double has only 53-bit resolution (e.g. 2^53+1 == 2^53)

numCol = ceil(N/numBits); % number of columns

stringsDecimal = zeros(2^20,1);
count = numCol; % include all zero states

for ind = 1:numel(maximalIS)
    mIS = maximalIS{ind}; % a list of node indices for mIS, e.g. [1,4,6,11]
    mISsize = ISsizes(ind);
    
    numNewIS = (2^mISsize-1);
    
    binaryFactor = 2.^(mIS(mIS<= numBits) - 1)';
    for ind2 = 2:numCol
        mIS = mIS(mIS>numBits);
        mIS = mIS - numBits;
        binaryFactor = blkdiag(binaryFactor, 2.^(mIS(mIS<= numBits) - 1)');
    end
    
 
    stringsNew = (binaryString{mISsize}*binaryFactor).';

    stringsDecimal(count+1:count+numCol*numNewIS) = stringsNew(:);
    
    count = count + numNewIS*numCol;
    
end

stringsDecimal(count+1:end) = [];
stringsDecimal = reshape(stringsDecimal, numCol,[])';


stringsDecimal = unique(stringsDecimal, 'rows');
numIS = size(stringsDecimal,1);

%%
stringsLegal = d2b(stringsDecimal.', numBits);
stringsLegal = reshape(stringsLegal.',[numCol*numBits, numIS]).';

stringsLegal = stringsLegal(:, 1:N);


function s=d2b(d,n)
% Adapted from DEC2BIN function

d = d(:); % Make sure d is a column vector.

% Actual algorithm
s = logical(rem(floor(d*pow2(0:-1:1-n)),2));
end


end

