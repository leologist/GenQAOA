function [NComp, wGraphComp, vertexComp] = CreateConnGraph(N,wGraph)
% This function divide the original graph into connected pieces
% Each with vertex labelled and the corresponding edges wGraphComp, and the number of vertex NComp

G = graph(wGraph(:,1),wGraph(:,2),[],N);
comp = conncomp(G);
NoComp = max(comp);

wGraphComp = cell(NoComp,1); % the wGraph in each connected piece, with vertex relabelled 
NComp = zeros(NoComp,1);
vertexComp =  cell(NoComp,1); % the vertices in the original label, in the corresponding order with wGraphComp

for indComp = 1:NoComp
    vertexCompT = find(comp == indComp); % connected vertex
    rowComp = ismember(wGraph(:,1), vertexCompT); % find the edges in the connected piece
    wGraphCompT = wGraph(rowComp,1:2); % edges in the connected piece
    [~,wGraphCompT] = ismember(wGraphCompT,vertexCompT); % reduce from the original vertex label to a relabelling in the connected piece, in the corresponding order with vertexCompT 
    wGraphComp{indComp} = wGraphCompT; 
    NComp(indComp) = numel(vertexCompT);
    vertexComp{indComp} = vertexCompT;
end