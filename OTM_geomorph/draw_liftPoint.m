%% load target mesh
model3 = createpde(1);
importGeometry(model3,'std_s_2rings.STL');
generateMesh(model3,'GeometricOrder','linear','Hmin',1,'Hmax',2);
pdeplot(model3);
n3=model3.Mesh.Nodes;
T3=model3.Mesh.Elements;
n3t=n3';
TR = triangulation(T3',n3t(:,1),n3t(:,2));
%%

model4 = createpde(1);
importGeometry(model4,'std_s_2circles.STL');
generateMesh(model4,'GeometricOrder','linear','Hmin',1,'Hmax',2);
pdeplot(model4);

n4=model4.Mesh.Nodes;
T4=model4.Mesh.Elements;
n4t=n4';
TS = triangulation(T4',n4t(:,1),n4t(:,2));      
%%
E=TS.Points;
wts=zeros(size(E,1),1);

disp('Lifting points');
LE = liftPD(E, wts);
disp('Computing convex hull');
C = convhulln(LE);
disp('Extracting lower hull');
ind = normalsPD(LE, C);
T = C(ind, :);

trisurf(T,LE(:,1),LE(:,2),LE(:,3));
hold on;
triplot(T,E(:,1),E(:,2));
% c_table=['k','b','g','y','r','m','c','w'];
% for k=1:1:test_cell_n
%     if ~isempty(PD{k, 1})        
%         tes=PD{k, 1};
%         tc=Polyhedron(tes);
%         col=mod(k,length(c_table))+1;
%         plot(tc,'color',c_table(col),'shade',0.1);
% %         plot(tc,'wire',1);
% %         kk=convhulln(tes);
% %         trisurf(kk,tes(:,1),tes(:,2),tes(:,3),'FaceColor','cyan','FaceAlpha',0.01)
% %         plot(tes(tc',1),tes(tc',2),'b-','LineWidth',2.6)
%         hold on;
%     end
% end