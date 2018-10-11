function M = krons(A1, A2, varargin)
%krons(A1, A2, ...) like kron() but accepts arbitrarily many arguments
%   returns A1 \otimes A2 \otimes ....
%   
%   M = kron(A1, A2, ...)
%
%   author: Leo Zhou
%   date: 4/17/2017

if ~isempty(varargin)
    M = krons(kron(A1,A2), varargin{1}, varargin{2:end});
else
    M = kron(A1,A2);
end

end

