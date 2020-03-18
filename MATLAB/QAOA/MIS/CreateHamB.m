function HamB = CreateHamB(stringsLegal, vecLegality, Omega)
%CreateHamB generates the Hamiltonian for \sum_i Omega_i*X_i in the
%           subspace of legal strings 
%
% Usage:
%   HamB = CreateHamB(stringsLegal, vecLegality)
%   HamB = CreateHamB(stringsLegal, vecLegality, Omega)
%
% Input:
%   stringsLegal, vecLegality = output of CreateSubspace function
%       stringsLegal = a KxN table indicating all independent sets (legal
%                     strings)
%
%       vecLegality = 2^N x 1 vector, where nonzero numbers occupy the 
%            independent set states in the original Hilbert space, each 
%            labelling the row position of the legal string in the order 
%            specified in the stringsLegal table.
%   
%   Omega = array of individual site coupling strength (default=ones(N,1))

NoString = size(stringsLegal,1);
N = size(stringsLegal,2);

if nargin == 2    
    Omega = ones(N,1); % homogeneous onsite strength
end

rowInd = zeros(NoString*N,1);
colInd = zeros(NoString*N,1);
HamBval = zeros(NoString*N,1);
binaryFactor = pow2(0:N-1)';

count = 1;
for indString = 1:NoString
    string = stringsLegal(indString,:);
    stringFlip = repmat(string,N,1); % repmat

    stringFlip(1:(N+1):end) = ~string; % single spin flip
    stringFlipind = stringFlip*binaryFactor + 1;
    [stringFlipSiteind,~,stringFlipind] = find(vecLegality(stringFlipind)); 
    % stringFlipSiteind: find the site index that has a connected string
    % stringFlipind: find the index in the reduced subspace
    stringFlipSize = numel(stringFlipind);
    rowInd(count:count+stringFlipSize-1) = indString;
    colInd(count:count+stringFlipSize-1) = stringFlipind;
    HamBval(count:count+stringFlipSize-1) = Omega(stringFlipSiteind);
    count = count + stringFlipSize;
end

rowInd(count:end) = []; colInd(count:end) = []; HamBval(count:end) = []; 
HamB = sparse(rowInd,colInd,HamBval,NoString,NoString);
