function H_AB = Ham2LTerm(sA, sB, iA, iB, n)
%Ham2LTerm(sA, sB, iA, iB, n) generate 2-local term in n-particle space
%   sA, sB: the two operators on qudit A and B
%   iA, iB: the indices of the two qudit
%   n: total number of qubits
%

if iA == iB
    error('iA cannot equal iB');
elseif iA < 0 || iB < 0 || n < max(iA, iB)
    error('invalid input of iA, iB, or n');
elseif ~isequal(size(sA), size(sB))
    error('operator sA and sB must be of same size');
end

d0 = length(sA);
if issparse(sA)
    myeye = @(n) speye(n);
elseif size(sA,2) == 1
    myeye = @(n) ones(n,1);
else
    myeye = @(n) eye(n);
end


if iA > iB
    H_AB = Ham2LTerm(sB, sA, iB, iA, n);
else
    H_AB = kron(sA, kron(myeye(d0^(iB-iA-1)), sB));
    H_AB = kron(myeye(d0^(iA-1)), kron(H_AB, myeye(d0^(n-iB))));
end

end

