function [PD,PDv] = pwrDiagramPD(T, PC, E)
% function PD = pwrDiagramPD(T, PC)
%
% T: triangulation
% PC: power centers of triangles in T
%
% The output cell PD contains pieces of the power diagram, indexed by
% dimension. PD{mP} contains pieces of dimension zero (power centers that
% correspond to fully-dimensional pieces of the triangulation T) and PD{1}
% contains fully-dimensional regions of the power diagram (corresponding to
% vertices of the triangulation).
P = piecesPD(T);
mP = size(P,1);
% count the empty cells for discrete points even they not occupy regions
PD = cell(mP,1);
PDv= cell(mP,1);
for i=1:mP
    EA = edgeAttPD(T, P{i}, E);
    for j=1:size(EA,1)
        PDv{i}{j,1}=EA{j};
        EA{j} = PC(EA{j},:);
    end
    PD{i} = EA;
end