function out = MultSingleSpin(psi,N,i,si)
% Action of I_2 tensor I_2 tensor ... s_i tensor I_2 ... I_2 on a state psi
% Do not form the tensor matrix explictly to greatly speed up the operation.
% See hTucker toolbox for a lot more tensor manipulations
% Only works for sx, sy, sz here for further speedup.
% For general si matrix, one can reshape and permute, back and forth

IndL = 2^(i-1); % i-1 is the number of sites on the left.
IndR = 2^(N-i); % on the right

% reshape for matrix operation
psi = reshape(psi, [IndR, 2, IndL]); %---- NOTE: the tensor product counts from right first

switch si
    case 1 % sx
        out = flip(psi,2);
    case 2 % sy
        out = zeros(size(psi));
        out(:,1,:) = -1i*psi(:,2,:);
        out(:,2,:) = 1i*psi(:,1,:);
    case 3 % sz
        out = psi;
        out(:,2,:) = -out(:,2,:);
end

out = reshape(out,[],1);

end