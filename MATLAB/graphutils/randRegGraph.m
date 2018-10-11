function [edges, A] = randRegGraph(n, d)
% RANDREGGRAPH(n,d) - creates a simple, approximately d-regular undirected graph
% simple = without loops or double edges
% d-reglar = each vertex is adjecent to d edges
% approximately d-regular: when n*d is not even, one vertex is allowed to
% have d-1 edges
%
% input arguments :
%   n - number of vertices
%   d - the degree of each vertex
% output :
%   edges - nx2 matrix specifying the edges
%   A - adjacency matrix of the generated graph

% check parameters
if n <= 0
    error(sprintf('invalid number of vertices n=%d',n));
end
if d >= n || d < 0
    error(sprintf('invalid degree d=%d. must be nonnegative and less than n',d));
end

if mod(n*d, 2) == 0
    edgesGoal = n*d/2;
else
    edgesGoal = (n*d-1)/2;
end

maxIter = 20; % how many attempts we want to try generating this graph

%a list of open half-edges
U = repmat(1:n,1,d);

%the graphs adajency matrix
A=sparse(n,n);

edgesTested=0; 
attempts=1;

%continue until a proper graph is formed
while ~isempty(U) && attempts <= maxIter

    %print progress
    if mod(edgesTested, 5000)==0  && edgesTested ~= 0
        currentEdges = nnz(A)/2;
        fprintf('randRegGraph() progress: attempt #%d/%d\n--- edges tried=%d, edges=%d/%d\n',...
            attempts, maxIter, edgesTested, currentEdges, edgesGoal);
        msg=isRegularGraph(A,d);
        
        if isempty(msg)
            break
        end
%         fprintf('Current graph%s\n', msg);
        
        %restart process if needed
        if (edgesTested > n*d*10 && edgesGoal-currentEdges <= 1) || (edgesTested > n*d*100) % have tested too many times
            attempts=attempts+1;            
            edgesTested = 0;
            U = repmat(1:n,1,d);
            A = sparse(n,n);
        end
    end

    %chose at random 2 half edges
    i1 = ceil(rand*length(U));
    i2 = ceil(rand*length(U));
    v1 = U(i1);
    v2 = U(i2);

    %check that there are no loops nor parallel edges
    if ~((v1 == v2) || (A(v1,v2) == 1))
        %add edge to graph
        A(v1, v2)=1;
        A(v2, v1)=1;
%         fprintf('found-%i,%i,%i-',v1,v2,full(A(v1,v2)));
        
        %remove used half-edges
        v = sort([i1,i2]);
        U = [U(1:v(1)-1), U(v(1)+1:v(2)-1), U(v(2)+1:end)];
    end
    
    edgesTested = edgesTested + 1;
end

% if ~isempty(U)
if length(U) > 1
    error(sprintf('Failure to generate a random regular graph after %d attempts', attempts));
end

[row, col] = find(triu(A));
edges = [row, col];
%-------------------------------------------------

function msg=isRegularGraph(G,d)
%is G a simple d-regular graph the function returns []
%otherwise it returns a message describing the problem in G

msg=[];

%check symmetry
if (norm(G-G','fro')>0)
    msg=[msg,' is not symmetric, '];
end

%check parallel edged
if (max(G(:))>1)
    msg=[msg,sprintf(' has %d parallel edges, ',length(find(G(:)>1)) )];
end

%check that d is d-regular
d_vec=sum(G);
if min(d_vec)<d-1 || max(d_vec)>d || length(find(d_vec == d-1)) > 1
    msg=[msg,' not approx d-regular, '];
end

%check that g doesn't contain any loops
if (norm(diag(G))>0)
    msg=[msg,sprintf(' has %d self loops, ',length(find(diag(G)>0)) )];
end