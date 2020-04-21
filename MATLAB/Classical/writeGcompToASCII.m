function writeGcompToASCII(filename, Adj)
%%writeGcompToASCII writes the complement graph specified by the given
%   adjacency matrix Adj to a text file in the DIMACS ASCII format
%

    fileID = fopen(filename, 'w');

    Nvert = length(Adj);
    Acomp = tril(~Adj,-1);
    Nedge = nnz(Acomp);
    
    fprintf(fileID, 'p col %d %d\n', Nvert, Nedge);
    for ind = 1:Nvert
        Js = find(Acomp(ind,:));
        if ~isempty(Js)
            for ind2 = Js
                fprintf(fileID, 'e %d %d\n', ind, ind2);
            end
        end
    end
    fclose(fileID);