function writeIsingToASCII(filename, h, J)
%%writeIsingToASCII writes Ising minimization problem to MaxSAT in an ASCII
% file in the DIMACS format
%
%   Ising problem: minimize sum_i h_i x_i + sum_{i<j} J_{ij} x_i x_j
%      where x_i \in {0,1}.
%
%   Equivalent weighted CNF for MaxSAT:
%       +x     => ~x (true if x = 0)
%       +x1 x2 => ~x1 | ~x2 (true iff x1 = x2 = 0)
%       -x1 x2 => (x1 | x2) & (x1 | ~x2) & (~x1 | x2) 
%                   ^ 3 clauses are true if x1 = x2 = 1, 2 otherwise
%       
%   
%   Usage:
%       writeIsingToASCII(filename, h, J)
%
%       h is length-N vector encoding h_i = h(i)
%       J is table, where each row is [i, j, J_{ij}]
%
%   DIMACS file format for MaxSat: http://www.maxsat.udl.cat/14/requirements/

fileID = fopen(filename, 'w');

numVars = length(h);
J_ij = J(:, 3);
numClauses = nnz(h) + nnz(J_ij > 0) + 3*nnz(J_ij < 0);

fprintf(fileID, 'p wcnf %d %d\n', numVars, numClauses);

% writing on-site field to clauses
for ind = 1:numVars
    if h(ind) ~= 0
        fprintf(fileID, '%d %d 0\n', abs(h(ind)), -sign(h(ind))*ind);
    end
end

% writing couplings to clauses

for row = 1:size(J,1)
    if J(row, 3) > 0
        fprintf(fileID, '%d %d %d 0\n', J(row, 3), -J(row,1), -J(row,2));
    elseif J(row, 3) < 0
        fprintf(fileID, '%d %d %d 0\n', -J(row, 3), J(row,1), J(row,2));
        fprintf(fileID, '%d %d %d 0\n', -J(row, 3), -J(row,1), J(row,2));
        fprintf(fileID, '%d %d %d 0\n', -J(row, 3), J(row,1), -J(row,2));
    end
end

fclose(fileID);

end