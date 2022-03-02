function [nPD, nPDv, nPD_vertex, cell_contri]=compute_restricted_polytopes(TS,PD)

cell_num=size(TS.ConnectivityList,1);
maxP_per_cell=size(TS.Points,1);
test_cell_n=size(PD{1,1},1);

nPD=cell(test_cell_n,1);


cell_contri=zeros(test_cell_n,1);
PDcell=PD{1, 1};

% it seems that the bound of destination mesh should be strictly convex
% give certain resolution but in current case even the set of cut cells is
% not convex the extracted extremes are used to comprise new large convex

parfor nn=1:1:test_cell_n 
    % set a number that would never be encountered within constrain mesh
    extreme_cell=zeros(maxP_per_cell+size(PDcell{nn, 1},1),2);
    ck=0;
for m=1:1:cell_num
        P1=PDcell{nn, 1};
        P2=TS.Points(TS.ConnectivityList(m,:),:);
        if (~isempty(P1))&&(~isempty(P2))
            I1=polytope(P1);
            I2=polytope(P2);
            % rough estimate of the cut-cell
            cut_cell=intersect(I1,I2);
            %pe=round(extreme(cut_cell),5);
            % the correct operation should be the truncating
            % in testing the Y aix has larger variance for cut (depending on specific geometry and slope of polytopes)
            % we can explore the different resolution for different stages
            % of optimization and scales of the geometry
            cc_p=extreme(cut_cell);
            pe=[fix(cc_p(:,1)*(10^5))/(10^5),fix(cc_p(:,2)*(10^4))/(10^4)];
            if (~isempty(pe))&&(ck==0)
            init_cpl=size(pe,1);
            extreme_cell(ck+1:ck+init_cpl,:)=pe;
            ck=ck+init_cpl;
            end

            if (~isempty(pe))&&(ck~=0)
            im=~ismember(pe,extreme_cell(1:ck,:),'rows');
            extreme_cell(ck+1:ck+sum(im),:)=pe(im,:);
            ck=ck+sum(im);
            end
            cell_contri(nn)=cell_contri(nn)+volume(cut_cell);

      
        end
end
I_cut=polytope(extreme_cell(1:ck,:));
nPD{nn}=extreme(I_cut);
end

sizeC=cellfun(@length, nPD);
nPD_vertex=unique(vertcat(nPD{:}),'rows');
[~,lb]=ismember(vertcat(nPD{:}),nPD_vertex,'rows');
nPDv=mat2cell(lb,sizeC);
end