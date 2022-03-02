clc;
clear all;
%% load 3d model from CAD
model = createpde(1);
importGeometry(model,'2sphere.stl');
generateMesh(model,'GeometricOrder','linear','Hmin',23);
pdeplot3D(model);

n1=model.Mesh.Nodes;
T1=model.Mesh.Elements;
T1t=T1';
n1t=n1'*0.05;
% n1t=n1';
tetramesh(T1t,n1t,'FaceAlpha',0.1);
TR=triangulation(T1t,n1t);
view(-40,20);
%% compute initial pd of the above with zero weights
E=TR.Points;
wts=zeros(size(E,1),1);
% wts=rand(size(E,1),1)*5;
[PD, PDinf, Votex, PDv]=PowerDiagramFunc(E, wts);
%% load src model
model2 = createpde(1);
importGeometry(model2,'1sphere.STL');
generateMesh(model2,'GeometricOrder','linear','Hmin',26);
pdeplot3D(model2);
n2=model2.Mesh.Nodes;
T2=model2.Mesh.Elements;
T2t=T2';
n2t=(n2'+[40,0,0])*0.05;
TS=triangulation(T2t,n2t);
%% compute volume of src
cell_num=size(TS.ConnectivityList,1);
area_cell=zeros(cell_num,1);
cell_cent=zeros(cell_num,2);
area_sum=0;
for i=1:1:cell_num
cell_tri=TS.Points(TS.ConnectivityList(i,:),:);
area_cell(i)=volume(Polyhedron(cell_tri));
cell_cent(i,:)=[mean(cell_tri(:,1)),mean(cell_tri(:,2))];
area_sum=area_sum+area_cell(i);
end
%% compute constrained voronoi diagram
tic
[nPD, nPDv, nPD_vertex, cell_contri]=compute_restricted_polytopes_3d(TS,PD);
toc
%% cut test
test_cell_n=size(E,1);
wts=zeros(size(E,1),1);
Nu=area_sum/test_cell_n;
Nu_p=zeros(test_cell_n,1);
Nu_p(:)=Nu;
Loop=500;
lr=0.9666666;
tStart=tic;
for i=1:1:Loop
    disp(['Loop',' ',num2str(i),' ','begins']);
    [PD, PDinf, vertex, PDv] = PowerDiagramFunc(E, wts);
    tst=tic;
    [nPD, nPDv, nPD_vertex, cell_contri]=compute_restricted_polytopes_3d(TS,PD);
    tend=toc(tst);
    gradient=-cell_contri+Nu_p;
    hessian=compute_hessian_3d(E, nPD, nPDv, nPD_vertex);
    disp('get direction');
    dh=-hessian\gradient;
    wts=wts+lr*dh;
    disp('Current norm of gradient is')
    norm(gradient)
    disp('One loop takes time of ')
    tend
end
tEnd=toc(tStart);
tEnd
%% do final cut and recompute restricted DT and get final connectivity matrix of Tet(G)
BD_p=TS.Points;
[Kp,total_mass]=convhulln(BD_p);
BD_pt=BD_p(Kp(:,1),:);
omitl=zeros(test_cell_n,1);
m_v_cent=zeros(test_cell_n,3);
for kp=1:1:test_cell_n
    try 
    m_v_cent(kp,:)=compute_centroid(nPD{kp});
    catch
    kp
    end
end
tablel=[1:test_cell_n]';
tablelt=[1:test_cell_n]';
tablel=tablel(omitl==0);
tablelt=tablelt(omitl==1);
V_init=m_v_cent(tablel,:);
v_init_bd=[m_v_cent;BD_pt];

% Tf=delaunayn(m_v_cent(tablel,:));
% Tg is the initial mesh triangulation
% Tg=T(sum(ismember(T,tablelt),2)==0,:);
LE = liftPD(E, wts);
C = convhulln(LE);
ind = normalsPD(LE, C);
T = C(ind, :);
Tg=T;
P_end=E+[10,0,0];
P_init=m_v_cent;

% the initial triangulation is customized in reality
T_end=TR.ConnectivityList;
Tet=intersect(sort(T_end,2),sort(Tg,2),'rows');
%% make GIF 
% h1 = figure(1);
loops = 26;
T_interval=1/loops;
% axis manual % this ensures that getframe() returns a consistent size

filename = 'newton3d_bm2_dynamic2_para5.gif';
for j = 1:1:(2*loops+1)
    h1 = figure(1);
    ax = axes();
    fig=gcf;
    fig.Position(3:4)=[1000,500];
%     get_3d_animation(T_end,P_t)
    xlabel('X');
    ylabel('Y');
    zlabel('Z');
    xlim(ax,[-5 20])
    ylim(ax,[-5 10])
    zlim(ax,[-5 10])

    if j<(1+loops)
    Tn=(j-1)*T_interval;
    P_t=Tn.*P_end+(1-Tn).*P_init;
    az=(j-1)*T_interval*175;
    view(ax, [az, 20])
    hold on;
    if j<3
        Tt=Tg;
    else
        Tt=T_end;
    end
    tetramesh(Tt,P_t);

    else
    Tn=(j-1-loops)*T_interval;
    P_t=Tn.*P_init+(1-Tn).*P_end;
    az=(j-1)*T_interval*175;
    view(ax, [az, 20])
    hold on;
    if j>2*loops-2
        Tt=Tg;
    else
        Tt=T_end;
    end
    tetramesh(Tt,P_t);

    end

    pbaspect([25 15 15]);
%     view(-45,20);
    
    drawnow
%     scatter(a2,P_t(:,1),P_t(:,2),'k*');
%     xlim([-10 120])
%     ylim([-10 120])
%     drawnow
      % Capture the plot as an image 
      frame = getframe(h1); 
      clf reset;
      im = frame2im(frame); 
      [imind,cm] = rgb2ind(im,256); 
      % Write to the GIF File 
      if j == 1 
          imwrite(imind,cm,filename,'gif','Loopcount',Inf); 
      else 
          imwrite(imind,cm,filename,'gif','WriteMode','append'); 
      end 
 end