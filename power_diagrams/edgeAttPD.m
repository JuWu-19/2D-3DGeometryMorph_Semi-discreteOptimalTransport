function EA = edgeAttPD(T, edges, E)
% function EA = edgeAttPD(T, edges)
%
% T: pieces of a triangulation
% edges: array of edges in the triangulation
% 
% Each entry in the output cell EA corresponds to an edge from the input
% edges and contains all pieces of the triangulation attached to that edge.
total_num=size(E,1);
[me, ne] = size(edges);
% EA = cell(me,1);
if ne==1
    EA = cell(total_num,1);
for i=1:1:total_num
    if ismember(i,edges)
    EA{i} = find(sum(ismember(T, i), 2) == ne)';
    else
    EA{i}=[];
    end
end
else
EA = cell(me,1);
for i=1:1:me
    EA{i} = find(sum(ismember(T, edges(i,:)), 2) == ne)';
end


end