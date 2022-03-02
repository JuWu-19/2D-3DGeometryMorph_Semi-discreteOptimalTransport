% sample input for 2D Voronoi diagram
clc;
clear all;
load('test_bm489p_110lps.mat');
%% make GIF 
h1 = figure(1);
loops = 18;
T_interval=1/loops;
axis tight manual % this ensures that getframe() returns a consistent size
filename = 'newton2e_bm1_rec2ring_dynamic1.gif';
for j = 1:1:(2*loops+1)

    fig=gcf;
    fig.Position(3:4)=[1000,300];
    if j<loops+1

    if j<2
        Tt=Tg;
    else
        Tt=T_end;
    end
    Tn=(j-1)*T_interval;
    P_t=Tn.*P_end+(1-Tn).*P_init;
    triplot(Tt,P_t(:,1),P_t(:,2))

    else

    if j>2*loops-1
        Tt=Tg;
    else
        Tt=T_end;
    end
    Tn=(j-1-loops)*T_interval;
    P_t=Tn.*P_init+(1-Tn).*P_end;
    triplot(Tt,P_t(:,1),P_t(:,2))




    end

    pbaspect([2 1 1])
    xlabel('X');
    ylabel('Y');
    xlim([0 200]);
    ylim([0 100]);
    drawnow
%     scatter(a2,P_t(:,1),P_t(:,2),'k*');
%     xlim([-10 120])
%     ylim([-10 120])
%     drawnow
      % Capture the plot as an image 
      frame = getframe(h1); 
      
      im = frame2im(frame); 
      [imind,cm] = rgb2ind(im,256); 
      % Write to the GIF File 
      if j == 1 
          imwrite(imind,cm,filename,'gif','Loopcount',Inf); 
      else 
          imwrite(imind,cm,filename,'gif','WriteMode','append'); 
      end 
 end