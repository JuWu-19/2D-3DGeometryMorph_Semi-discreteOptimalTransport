function H=compute_hessian_3d(E, PD, PDv, vertex)
test_cell_n=size(PD,1);
PD_neigh=cell(test_cell_n,1);
% get neighbors of each cell
for i=1:1:size(test_cell_n,1)
    PD_neigh{i}=[];
end

for k=1:1:test_cell_n-1
    comp_ini=k+1;
for j=comp_ini:1:test_cell_n
        if sum(ismember(PDv{k},PDv{j}))~=0
            if sum(ismember(PDv{k},PDv{j}))<3
                % this error indicator is of contradiction with rounds
                % constrained voronoi cells
                % it is hard to find out an universal cut resolution so it
                % maybe reasonable to dismiss pairs of 1 common point even
                % it would dismiss false negative ones (negative means only one common point between cells)
                disp('Common points error detected, they are points of cells');
                k
                j
                disp('They have the following number of points intersected')
                sum(ismember(PDv{k},PDv{j}))
                continue;
            end
            PD_neigh{k}=[PD_neigh{k},j];
            PD_neigh{j}=[PD_neigh{j},k];
        end

end 

end
%total number of facets
num_facet=sum(cellfun(@length,PD_neigh));
I=zeros(num_facet,1);
J=zeros(num_facet,1);
Vij=zeros(num_facet,1);
k=0;
for m=1:1:size(PDv,1)
    if isempty(PDv{m,1})
        continue;
    else
    for n=1:1:length(PD_neigh{m, 1})
    k=k+1;
    I(k)=m;
    J(k)=PD_neigh{m, 1}(n);
    neighbor_n=PDv{PD_neigh{m, 1}(n), 1};
    if size(neighbor_n,1)<4
        size(neighbor_n,1)
        disp('Mesh cut error detected');
    end
    cut_facet=intersect(neighbor_n,PDv{m,1});
    if length(cut_facet)<3
        disp('Only cut line error');
        area_cut=0;
    elseif length(cut_facet)==3
        area_cut=compute_surface_3d(vertex(cut_facet,:));
    else 
        try
            c_facet=convhulln(vertex(cut_facet,:));
            disp('Not coplanarity error');
            PD_neigh{m, 1}(n)
            m
            area_cut=compute_surface_3d(vertex(cut_facet,:));
        catch 
            disp('Multiple points coplanarity');
            area_cut=compute_surface_3d(vertex(cut_facet,:));
        end

    end

    Vij(k)=area_cut/pdist(E([m,PD_neigh{m, 1}(n)],:),'euclidean');
    end
    end
end
disp('compute hessian');
delta=-0.665;
Pert=delta*ones(test_cell_n,1);
Hp=sparse(diag(Pert));
size(Hp)
H1=sparse(I,J,Vij,test_cell_n,test_cell_n);
size(H1)
Vii=-sum(H1,2);
H2=sparse(diag(Vii));
H=H1+H2+Hp;
end