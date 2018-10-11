function K = krondist(M, k, varargin)
%KRONDIST returns an operator based on M distributed over k subsystems
%
%   K = KRONDIST(M, k) returns (MIII..+IMII..+IIMI+..), where I is the
%   identity matrix with the same dimension as M
%   For example, when k = 3, KRONDIST returns (MII+IMI+IIM), where
%   the products are implicitly kronecker products
%
%   K = KRONDIST(M, k, coefs) returns  $K=sum_{i=1}^k M(i)*coefs(i)$

if issparse(M)
    myeye = @(n) speye(n);
elseif size(M,2) == 1
    myeye = @(n) ones(n,1);
else
    myeye = @(n) eye(n);
end
d = length(M);

coefs = ones(1,k);

if ~isempty(varargin)
    if length(varargin{1}) ~= k
        error('length of coefficent vector inconsistent with k');
    end
    coefs = varargin{1};
end
    

K = sparse(0);
for ind = 1:k
    K = K + kron(kron(myeye(d^(ind-1)), M), myeye(d^(k-ind))) * coefs(ind);
end


end

