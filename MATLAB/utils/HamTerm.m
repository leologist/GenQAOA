function h = HamTerm(ops, sites, N)
%HamTerm(ops, inds, N)
%
%   Example:
%       h = HamTerms({X, Y, X}, [2, 3, 5], 6);
%     returns IXYIXI

if ~issorted(sites)
    [s, Isort] = sort(sites);
    h = HamTerm(ops(Isort), s, N);
    return
end

d0 = length(ops{1});
if issparse(ops{1})
    myeye = @(n) speye(n);
elseif size(ops{1},2) == 1
    myeye = @(n) ones(n,1);
else
    myeye = @(n) eye(n);
end


prev_site = 0;
h = 1;
for ind = 1:length(sites)
    site = sites(ind);
    h = krons(h,  myeye(d0^(site-prev_site-1)), ops{ind});
    prev_site = site;
end
h = kron(h, myeye(d0^(N-sites(end))));

end