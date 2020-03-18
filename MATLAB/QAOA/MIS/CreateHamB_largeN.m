function HamB = CreateHamB_largeN(stringsLegal, Omega)
%CreateHamB_largeN generates the Hamiltonian for \sum_i Omega_i*X_i in the
%           subspace of legal strings (independent sets), should be used
%           when N is large (N >= 48)
%
%   Usage: HamB = CreateHamB_largeN(stringsLegal)
%        HamB = CreateHamB_largeN(stringsLegal, Omega)
%
%   Input:
%        stringsLegal = a logical array of N columns representing
%                       bitstrings corresponding to independent sets
%
%        Omega(i) = Omega_i is the i-th site coupling strength
%            will assumed to be 1 if this argument is not supplied


NumString = size(stringsLegal,1);
N = size(stringsLegal,2);

if nargin < 2    
    Omega = ones(N,1); % homogeneous onsite strength
end

%% constructing mixing Hamiltonian B

rowInd = zeros(NumString*N,1);
colInd = zeros(NumString*N,1);
HamBval = zeros(NumString*N,1);

count = 1;

strHashes = hash_bitstrings_to_ints(stringsLegal);

for indString = 1:NumString
    string = stringsLegal(indString,:);
    stringFlipped = repmat(string,N,1); % pre-generate array for all possible spin flips

    stringFlipped(1:(N+1):end) = ~string; % perform single spin flip

    % logical index array of legal strings connected via single-site flips
    stringFlippedHashes = hash_bitstrings_to_ints(stringFlipped);
    IstringFlippedLegal = ismember(strHashes, stringFlippedHashes, 'rows');

    % Alternative method without hashing - slower:
    % IstringFlippedLegal = ismember(stringsLegal, stringFlipped, 'rows');

    % number of legal strings connected via single-site flips
    stringFlippedSize = nnz(IstringFlippedLegal);
    
    % the site indices of the flipped spin for each flipped string
    stringFlippedSiteInds = mod(find(xor(stringsLegal(IstringFlippedLegal, :), string).'), N);
    stringFlippedSiteInds(stringFlippedSiteInds == 0) = N;
    
    rowInd(count:count+stringFlippedSize-1) = indString;
    colInd(count:count+stringFlippedSize-1) = find(IstringFlippedLegal);
    HamBval(count:count+stringFlippedSize-1) = Omega(stringFlippedSiteInds);
    count = count + stringFlippedSize;
end

rowInd(count:end) = []; colInd(count:end) = []; HamBval(count:end) = []; 
HamB = sparse(rowInd,colInd,HamBval,NumString,NumString);

%%
    function integers = hash_bitstrings_to_ints(bitstr)
    %% This function convert rows of bitstrings into rows of integers,
    %  uniquely representing (hashing) each bitstring for more efficient
    %  comparison.

        maxNumBits = 50;

        totalNumBits = size(bitstr, 2);
        numCol = ceil(totalNumBits/maxNumBits);

        integers = zeros(size(bitstr,1), numCol);
        binaryFactor = 2.^(0:maxNumBits-1)';

        for ind = 1:numCol-1
            % convert bits to integer
            integers(:,ind) = bitstr(:, maxNumBits*(ind-1) + (1:maxNumBits))*binaryFactor;
            totalNumBits = totalNumBits - maxNumBits;
        end

        binaryFactor = 2.^(0:totalNumBits-1)';
        integers(:, numCol) = bitstr(:, maxNumBits*(numCol-1) + (1:totalNumBits))*binaryFactor;

    end

end

