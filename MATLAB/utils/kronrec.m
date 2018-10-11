function K = kronrec(M, k)
%KRONREC recursively generate a matrix that is the kronecker product
%   of the matrix M with itself k times. In other words,
%   KRONREC returns M^{\otimes k}
%
%   K = kronrec(M, k)

if (mod(k, 1) ~= 0) || (k < 0)
    error('KRONREC: Invalid input k = %d', k);
end

if k == 0
    K = 1;
elseif k == 1
    K = M;
else
    K = kron(M, kronrec(M, k-1));
end

end

