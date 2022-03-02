function C = compute_centroid(P)
[~,k]=convhulln(P);
if abs(k-volume(Polyhedron(P)))>1e-3
    length(unique(k(:)))
    size(P,1)
    error('Polyhedron is not convex.');
end
T = delaunayn(P);
n = size(T,1);
W = zeros(n,1);
C=0;
for m = 1:n
    sp = P(T(m,:),:);
    [null,W(m)]=convhulln(sp);
    C = C + W(m) * mean(sp);
end
C=C./sum(W);
return