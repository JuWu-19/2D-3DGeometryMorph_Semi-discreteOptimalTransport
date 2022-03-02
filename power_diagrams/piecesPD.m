function [P, total] = piecesPD(T)
% function [P, total] = piecesPD(T)
%
% T: triangulation
%
% The output P contains all pieces of the triangulation, from vertices
% (dimension zero) to fully-dimensional facets (dimension n-1).

[m, n] = size(T);
P = cell(n,1);
total = 0;

P{1} = unique(T);
P{1} = P{1}(:);
total = total + size(P{1},1);

for i=2:n-1
%     Q=[]
        Q = zeros(nchoosek(n,i)*size(T,1),i);
    for j=1:m
        ck=(j-1)*nchoosek(n,i);
        ce=j*nchoosek(n,i);
%         Q = [Q; combnk(T(j,:),i)];
        Q((ck+1):ce,:)=nchoosek(T(j,:),i);
    end
    Q = sort(Q, 2);
    P{i} = unique(Q, 'rows');
    total = total + size(P{i},1);
    clear Q;
end

P{n} = T;
total = total + size(P{n},1);