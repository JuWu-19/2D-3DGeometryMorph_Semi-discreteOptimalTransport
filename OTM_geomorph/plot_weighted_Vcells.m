function plot_weighted_Vcells(PD, PDv, vertex, E, wts)

figure(1)
% LE = liftPD(E, wts);
% C = convhulln(LE);
% ind = normalsPD(LE, C);
% T = C(ind, :);
% hold on;
% [PC, ~] = powercentersPD(T, E, wts);
% scatter(E(:,1),E(:,2),150,'g*');
% a = [1:size(E,1)]'; b = num2str(a); c = cellstr(b);
% dx = 0.1; dy = 0.1; % displacement so the text does not overlay the data points
% text(E(:,1)+dx, E(:,2)+dy, c);
% hold on;
% 
% a2 = [1:size(T,1)]'; b2 = num2str(a2); c2 = cellstr(b2);
% dx = 0.1; dy = 0.1; % displacement so the text does not overlay the data points
% text(PC(:,1)+dx, PC(:,2)+dy, c2);
% scatter(PC(:,1),PC(:,2), 260,'k*');

test_cell_n=size(PD,1);

if size(E,2)==2

a = [1:size(vertex,1)]'; b = num2str(a); c = cellstr(b);
dx = 0.1; dy = 0.1; % displacement so the text does not overlay the data points
text(vertex(:,1)+dx, vertex(:,2)+dy, c);
hold on;
scatter(vertex(:,1),vertex(:,2),150,'g*');
hold on;

for k=1:1:test_cell_n
    if ~isempty(PD{k, 1})
        tes=PD{k, 1};
        tc=convhulln(tes);
        plot(tes(tc',1),tes(tc',2),'b-','LineWidth',2.6)
        hold on;
    end
end

for k=1:1:test_cell_n
        if ~isempty(PDv{k, 1})
            sq=PDv{k, 1};
            tes1=vertex(sq,:);
            tc1=convhulln(tes1);
            plot(tes1(tc1',1),tes1(tc1',2),'r-','LineWidth',1.5)
            hold on;
        end
end

end

if size(E,2)==3
% LE = liftPD(E, wts);
% C = convhulln(LE);
% ind = normalsPD(LE, C);
% T=C(ind, :);
% 
% num_pc=size(T,1);
% a = [1:size(vertex,1)]'; b = num2str(a); c = cellstr(b);
% dx = 0.01; dy = 0.01; dz=0.01; % displacement so the text does not overlay the data points
% text(vertex(:,1)+dx, vertex(:,2)+dy, vertex(:,3)+dz, c);
% hold on;
% scatter3(vertex(1:num_pc,1),vertex(1:num_pc,2),vertex(1:num_pc,3),150,'g*');
% hold on;
% scatter3(vertex(num_pc+1:end,1),vertex(num_pc+1:end,2),vertex(num_pc+1:end,3),350,'r*');
% hold on;
c_table=['k','b','g','y','r','m','c','w'];
for k=1:1:test_cell_n
    if ~isempty(PD{k, 1})        
        tes=PD{k, 1};
        tc=Polyhedron(tes);
        col=mod(k,length(c_table))+1;
        plot(tc,'color',c_table(col),'shade',0.1);
%         plot(tc,'wire',1);
%         kk=convhulln(tes);
%         trisurf(kk,tes(:,1),tes(:,2),tes(:,3),'FaceColor','cyan','FaceAlpha',0.01)
%         plot(tes(tc',1),tes(tc',2),'b-','LineWidth',2.6)
        hold on;
    end
end

% for k=1:1:test_cell_n
%         if ~isempty(PDv{k, 1})
%             sq=PDv{k, 1};
%             tes1=vertex(sq,:);
%             tc1=convhulln(tes1);
%             plot(tes1(tc1',1),tes1(tc1',2),'r-','LineWidth',1.5)
%             hold on;
%         end
% end


end


end