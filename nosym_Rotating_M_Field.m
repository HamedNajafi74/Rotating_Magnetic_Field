 %% Clearing variables and close windoes and clearing command window   
 
    close all
    clear
    
    clc
    
%     tic;
%% File Name

    filename = 'NO_Ideal_Animated5.gif';       % Specify the output file name
    gifc = 0 ;                                   %1==yes or 0==no ==> creating .gif animation file
%% variables definition:

    j=1i;                      %j=[0+1i]
    f=50;                      %frequency
    P=1/f;                     %Period
    n=70;                      %Number of time steps
    t=linspace(0,2*P,n);       %time
    w=2*pi*f;                  %angular frequency(Omega)
    Bm=1;                      %amplitude of Mfield produced by each phase
    m=round(sqrt(n))+1;
%% Here we can make fault:

    R =  1  +  0.10*sin(w*t)   -  0.05 ;                  %Br=R*Bm*sin(w*t); 
    S =  1  +  0.19*sin(3*w*t)   -  1.00 ;                  %Bs=S*Bm*sin(w*t-2*pi/3);
    T =  1  +  0.10*sin(2*w*t)   -  0.20 ;                  %Bt=T*Bm*sin(w*t+2*pi/3);
   
                        A=1.0 ; B=1.0 ;                             %Bnet :: [theta,Bn]=cart2pol(A*Nx,B*Ny);
   
%     R=1;S=1;T=1;                             %
    
%% Magnetic Fields definition:

    Br=(Bm*R).*sin(w*t)          ;% +     0.10*(Bm*R).*sin(3*w*t);             %Mfield produced by   phase  'R'  color=Black  +   3rd Harmonic
    Bs=(Bm*S).*sin(w*t-2*pi/3)   ;% +     0.15*(Bm*S).*sin(3*(w*t-2*pi/3));    %Mfield produced by   phase  'S'  color=Blue   +   3rd Harmonic
    Bt=(Bm*T).*sin(w*t+2*pi/3)   ;% +     0.12*(Bm*T).*sin(3*(w*t+2*pi/3));    %Mfield produced by   phase  'T'  color=Pink   +   3rd Harmonic
    
%% Mfields in cartesian coordinates:

    [Rx,Ry]=pol2cart(0,Br);
    %Rx=Br;Ry=Br-Br;
    [Sx,Sy]=pol2cart(-2*pi/3,Bs);
    %Sx=Bs*cos(-2*pi/3);Sy=Bs*sin(-2*pi/3);
    [Tx,Ty]=pol2cart(2*pi/3,Bt);
    %Tx=Bt*cos(2*pi/3);Ty=Bt*sin(2*pi/3);
    
%% Avaluating Bnet vector:

    Nx=A.*(Rx+Sx+Tx);
    Ny=B.*(Ry+Sy+Ty);
    
%%  Bnet in polar coordinates

    [theta,Bn]=cart2pol(Nx,Ny);      %Bnet Mfield      color=Red
    
%% Wave shapes and geometrical positions

    c1=sin(w*t)+j*cos(w*t);
    c2=1.5*c1;
%     cr=max(R)*c1;
%     cs=max(S)*c1;
%     ct=max(T)*c1;
    cN=Nx+j*Ny;
    ir=(Bm).*sin(w*t);         
    is=(Bm).*sin(w*t-2*pi/3);   
    it=(Bm).*sin(w*t+2*pi/3);
    in=t-(-1.5+t);
    iz=in-1.5;
%% Somthing for creating Animation
    if (gifc == 1) 
        
        
        im=cell(m);
        
        frame=struct('cdata',[],'colormap',[]);
        
    end

%% 
       
%   rd=zeros(1,n);
%   sd=rd;
%   td=rd;
 %   tt=t-(t-1.5);
    
%% Openning a Figure   

    figure('NumberTitle','off','Name','Magnetic Fields in the Stator By H.Najafi','Units','centimeters','Position',[0.25 2 40 18]);
    
    %% Animated Lines
    
    subplot(3,4,[3,4]); 
    %plot([t,iz],'Color','k','LineWidth',0.5,[t,ir],'Color','g','LineWidth',1);
    plot(t,iz,'k',t,ir,'g');
    title('~Current of each phase (R,S,T) per unit');
    ral=animatedline('Color','k','LineWidth',1.5,'MaximumNumPoints',n);
    ylim([-1.5 1.5]);xlim([0 max(t)]);grid on;grid minor;
    ax = gca;
    ax.YLabel.String = 'R';
    ax.YLabel.FontSize = 10;
    
    
    subplot(3,4,[7,8]);
    plot(t,iz,'k',t,is,'g');
    sal=animatedline('Color','b','LineWidth',1.5,'MaximumNumPoints',n);
    ylim([-1.5 1.5]);xlim([0 max(t)]);grid on;grid minor;
    ax = gca;
    ax.YLabel.String = 'S';
    ax.YLabel.FontSize = 10;
    
    
    subplot(3,4,[11,12]);
    plot(t,iz,'k',t,it,'g');
    tal=animatedline('Color','m','LineWidth',1.5,'MaximumNumPoints',n);
    ylim([-1.5 1.5]);xlim([0 max(t)]);grid on;grid minor;
    ax = gca;
    ax.YLabel.String = 'T';
    ax.YLabel.FontSize = 10;
    ax.XLabel.String = 'Time';
    ax.XLabel.FontSize = 10;
    ax.XLabel.FontWeight='bold';
    
    subplot(3,4,[9,10]);
    plot(360*f*t,in,'g','LineWidth',0.75);
    nal=animatedline('Color','r','LineWidth',1.5,'MaximumNumPoints',n);
    %ylim([-2 2]);
    xlim([0 (180/pi)*w*max(t)]);grid on; grid minor;
    ax = gca;
    ax.YLabel.String = 'Amplitude of Bnet';
    ax.YLabel.FontSize = 10;
    ax.YLabel.FontWeight='bold'; 
    ax.XTickMode = 'manual';
    ax.XTick = linspace(0,360*f*max(t),13);
    ax.XLabel.String = 'Position of Bnet(0-360 Degree)';
    ax.XLabel.FontSize = 10;
    ax.XLabel.FontWeight='bold';
    
    
    subplot(3,4,[1,2,5,6]); 
    plot(cN,'r','LineWidth',2);
    title('Magnetic Fields');
    ax = gca;
    ax.XLabel.String = 'Hamed Najafipour (S.R.T.T.U.-1396)';ax.XLabel.FontSize = 10;ax.XLabel.FontWeight='bold';
    %axis square;
    axis([-1.7 1.7 -1.7 1.7]);
    axis equal;
    grid on;grid minor;hold on;
%%
    
for ii=1:n  
        
   %   rd(ii)=Br(ii);
   %   sd(ii)=Bs(ii);
   %   td(ii)=Bt(ii);
      
      subplot(3,4,[3,4]);
      
      %title("B produced by each phase");
      
      
      %cla;
   % plot(t,ir,'g','LineWidth',1);
     
     addpoints(ral,t(ii),Br(ii));
     drawnow ;
      
    %  hold on;
%       plot(t,ir,'g','LineWidth',2);
%       hold on;
%       plot(t,rd,'k','LineWidth',2);
  %    grid on;
     % hold on;
     
     
  %    ylim([-1.5 1.5]);
    %  hold on;
     
     % hold  off;
      
      subplot(3,4,[7,8]);
      
     addpoints(sal,t(ii),Bs(ii));
     drawnow ;
      
      %cla;
      
%       plot(t,is,'g',t,sd,'b','LineWidth',2);
%       hold on;
    %  plot(t,is,'g','LineWidth',2);
     % hold on;
     % plot(t,sd,'b','LineWidth',2);
     % grid on;
      
      
  %    ylim([-1.5 1.5]);
     % hold on;
     
     subplot(3,4,[11,12]);
     
     addpoints(tal,t(ii),Bt(ii));
     drawnow ;
      
      
      
  %    cla;
      
   %   plot(t,it,'g',t,td,'b','LineWidth',2);
   %   hold on;
%       plot(t,it,'g','LineWidth',2);
%       hold on;
%       plot(t,td,'m','LineWidth',2);
   %   grid on;
      
      
   %   ylim([-1.5 1.5]);
      %hold on;
      
       subplot(3,4,[9,10]);
      
      addpoints(nal,360*f*t(ii),Bn(ii));
      drawnow ;
      
        
     subplot(3,4,[1,2,5,6]);
      
     cla;
     
     plot(c1,'g','LineWidth',1);
     hold on;
     plot(c2,'g','LineWidth',1);
     hold on; 
     plot(cN,'r','LineWidth',1); 
      hold on;
     
      
      
      
%       title('Rotating M.Field By Hamed Najafi');
      
      %legend('R','S','T','Bnet');
      
      plot([0 Nx(ii)],[0 Ny(ii)],'r','Marker','s','LineWidth',1.5);
      hold on;
      plot([0 Rx(ii)],[0 Ry(ii)],'k','Marker','s','LineWidth',1);
      hold on;
      plot([0 Sx(ii)],[0 Sy(ii)],'b','Marker','s','LineWidth',1);
      hold on;
      plot([0 Tx(ii)],[0 Ty(ii)],'m','Marker','s','LineWidth',1);
      hold on;
      
      
%       grid on;
     
%       axis square;
%       axis([-1.7 1.7 -1.7 1.7]);
      
      %legend('r=1','r=1.5','R','S','T','Bnet');
      
     % drawnow
     if (gifc == 1) 
         
                 frame(ii) = getframe(1);
                   im{ii} = frame2im(frame(ii));
     end
      % hold off;
      
     
      
end
%% Movie
%     close;
%     fig=figure('NumberTitle','off','Name',['Magnetic Fields in the Stator',date],'Units','centimeters','Position',[0.25 2 40 18]);
%     movie(fig,frame,2);
%% These following codes are written for exporting plots to a .gif (animation) file

  if (gifc == 1) 
        %close;
        figure;
        nImages=length(t);
        for idx = 1:nImages
            subplot(m,m,idx)
            imshow(im{idx});
        end


        for idx = 1:nImages
            [A,map] = rgb2ind(im{idx},256);
            if idx == 1
                imwrite(A,map,filename,'gif','LoopCount',Inf,'DelayTime',0);
            else
                imwrite(A,map,filename,'gif','WriteMode','append','DelayTime',0);
            end
        end
        close all;
  end

%%
%     toc;
%% The End... ...!?
