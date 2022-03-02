function [nPD, nPDv, nPD_vertex, cell_contri]=compute_restricted_polytopes_3d(TS,PD)

cell_num=size(TS.ConnectivityList,1);
maxP_per_cell=size(TS.Points,1);
test_cell_n=size(PD{1,1},1);

nPD=cell(test_cell_n,1);


cell_contri=zeros(test_cell_n,1);
PDcell=PD{1, 1};

PDcesll_H=cell(test_cell_n,1);
PDcell_H_ind=zeros(test_cell_n,1);
TS_H=cell(cell_num,1);
TS_H_ind=zeros(cell_num,1);
parfor nn=1:1:test_cell_n 
    % set a number that would never be encountered within constrain mesh
    P1=PDcell{nn, 1};
    if ~isempty(P1)
    PDcell_H{nn}=Polyhedron(P1).computeHRep;
    else
    PDcell_H_ind(nn)=1;
    end
end

parfor m=1:1:cell_num
    % set a number that would never be encountered within constrain mesh
    P2=TS.Points(TS.ConnectivityList(m,:),:);
    if ~isempty(P2)
    TS_H{m}=Polyhedron(P2).computeHRep;
    else
    TS_H_ind(m)=1;
    end
end
% it seems that the bound of destination mesh should be strictly convex
% give certain resolution but in current case even the set of cut cells is
% not convex the extracted extremes are used to comprise new large convex
% in the future
parfor nn=1:1:test_cell_n 
    % set a number that would never be encountered within constrain mesh
    extreme_cell=zeros(maxP_per_cell+size(PDcell{nn, 1},1),3);
    ck=0;
for m=1:1:cell_num
%         P1=PDcell{nn, 1};
%         P2=TS.Points(TS.ConnectivityList(m,:),:);
%         if (~isempty(P1))&&(~isempty(P2))
%             I1=Polyhedron(P1);
%             I2=Polyhedron(P2);
%             cut_cell=intersect(I1,I2);

        if (TS_H_ind(m)==0)&&(PDcell_H_ind(nn)==0)
            cut_cell=intersect(PDcell_H{nn},TS_H{m});

%             try
%             I1=Polyhedron(P1);
%             I2=Polyhedron(P2);
%             catch
%             P1
%             P2
%             end
%             % rough estimate of the cut-cell
%             try
%             cut_cell=intersect(I1,I2);
%             catch
%             continue;
%             end
            %pe=round(extreme(cut_cell),5);
            % the correct operation should be the truncating
            % in testing the Y aix has larger variance for cut (depending on specific geometry and slope of polytopes)
            % we can explore the different resolution for different stages
            % of optimization and scales of the geometry
            cc_p=cut_cell.V;
            pe=[fix(cc_p(:,1)*(10^5))/(10^5),fix(cc_p(:,2)*(10^5))/(10^5),fix(cc_p(:,3)*(10^5))/(10^5)];
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
            if size(cut_cell.V,1)>3
                try
                    m_cell=volume(cut_cell);
                    if isinf(m_cell)||isnan(m_cell)
                        disp('Inf volume error')
                        continue;
                    end
                    cell_contri(nn)=cell_contri(nn)+m_cell;
                catch
                    continue;
                end
            end

      
        end
%         clear P1 P2 I1 I2 cut_cell;
end
I_cut=Polyhedron(extreme_cell(1:ck,:));
nPD{nn}=I_cut.V;
%     clear extreme_cell;
end

sizeC=cellfun(@length, nPD);
nPD_vertex=unique(vertcat(nPD{:}),'rows');
[~,lb]=ismember(vertcat(nPD{:}),nPD_vertex,'rows');
nPDv=mat2cell(lb,sizeC);
end