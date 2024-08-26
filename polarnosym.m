%% Clearing variables and close windoes and clearing command window

close all
clear
clc

tic;
%% 

%opengl hardware;

%% File Name

filename = 'polarnosym41.gif';       % Specify the output file name

gifc = 0 ;                                    % 1==yes or 0==no ==> creating .gif animation file

%% variables definition:

j=1i;                      %j=[0+1i]
f=50;                      %frequency
P=1/f;                     %Period
n=100;                      %Number of time steps
t=linspace(0,2*P,n);       %time
    %w=2*pi*f                  %angular frequency(Omega) (rad/s)
Bm=1;                      %amplitude of Mfield produced by each phase
m=round(sqrt(n))+1;
wtd=360*f*t;               %[omega]*t (degree)
wtr=2*pi*f*t;              %[omega]*t (radian)

%% Here we can make fault:

R =  1  +  0.10*sin(1.*wtr)   -  0.05 ;                  %Br=R*Bm*sin(w*t);
S =  1  +  0.19*sin(3.*wtr)   -  0.20 ;                  %Bs=S*Bm*sin(w*t-2*pi/3);
T =  1  +  0.10*sin(2.*wtr)   -  0.20 ;                  %Bt=T*Bm*sin(w*t+2*pi/3);
                        %%A=1.0 ; B=1.0 ;                             %Bnet :: [theta,Bn]=cart2pol(A*Nx,B*Ny);
  %   R=1;S=1;T=1;                             %

%% Magnetic Fields definition:

Br=(Bm*R).*sin(wtr)          ;% +     0.10*(Bm*R).*sin(3*wtr);             %Mfield produced by   phase  'R'  color=Black  +   3rd Harmonic
Bs=(Bm*S).*sin(wtr-2*pi/3)   ;% +     0.15*(Bm*S).*sin(3*(wtr-2*pi/3));    %Mfield produced by   phase  'S'  color=Blue   +   3rd Harmonic
Bt=(Bm*T).*sin(wtr+2*pi/3)   ;% +     0.12*(Bm*T).*sin(3*(wtr+2*pi/3));    %Mfield produced by   phase  'T'  color=Pink   +   3rd Harmonic

%% Mfields in polar coordinates:

RP=complex(Br,0);           
SP=Bs.*exp(j.*(-2*pi/3));              
TP=Bt.*exp(j.*(2*pi/3));           

%% Avaluating Bnet(PBn) vector:

PBn=RP+SP+TP;
theta=angle(PBn);
Bn=abs(PBn);

%% Wave shapes and geometrical positions

c1=Bm.*exp(j.*wtr);
c2=1.5*c1;
ir=(Bm).*sin(wtr);              %Ideal wave form of R
is=(Bm).*sin(wtr-2*pi/3);       %Ideal wave form of S
it=(Bm).*sin(wtr+2*pi/3);       %Ideal wave form of T
in=t-(-1.5+t);
iz=zeros(1,n);

%% Pre-allocating MField Vectors
a(n)=line;                  %a,b,c are M.field Vectors
for ii=1:n      
    a(ii)=line;  
end
b=a;
c=b;
d=c;
close;
%% Somthing for creating Animation

if (gifc == 1)    
    im=cell(n,1);    
    frame=struct('cdata',n,'colormap',3);
end

%%

%   rd=zeros(1,n);
%   sd=rd;
%   td=rd;
%   tt=t-(t-1.5);

%% Openning a Figure

figure('Renderer','OpenGL','NumberTitle','off','Name','Magnetic Fields in a 3phase-Stator By H.Najafi','Units','centimeters','Position',[0.25 2 40 18]);

%% Animated Lines

subplot(3,4,[3,4]);
plot(t,iz,'k',t,ir,'g');
title('~Current of each phase (R,S,T) per unit');
ral=animatedline('Color','k','LineWidth',1.5,'MaximumNumPoints',n);
ylim([-1.5 1.5]);xlim([0 max(t)]);
grid on;
grid minor;
ax = gca;
ax.YLabel.String = 'R';
ax.YLabel.FontSize = 10;
ax.YLabel.FontWeight='bold';

subplot(3,4,[7,8]);
plot(t,iz,'k',t,is,'g');
sal=animatedline('Color','b','LineWidth',1.5,'MaximumNumPoints',n);
ylim([-1.5 1.5]);xlim([0 max(t)]);
grid on;
grid minor;
ax = gca;
ax.YLabel.String = 'S';
ax.YLabel.FontSize = 10;
ax.YLabel.FontWeight='bold';

subplot(3,4,[11,12]);
plot(t,iz,'k',t,it,'g');
tal=animatedline('Color','m','LineWidth',1.5,'MaximumNumPoints',n);
ylim([-1.5 1.5]);xlim([0 max(t)]);
grid on;
grid minor;
ax = gca;
ax.YLabel.String = 'T';
ax.YLabel.FontSize = 10;
ax.YLabel.FontWeight='bold';
ax.XLabel.String = 'Time';
ax.XLabel.FontSize = 10;
ax.XLabel.FontWeight='bold';

subplot(3,4,[9,10]);
plot(wtd,in,'g','LineWidth',0.75);
title('Hamed Najafipour (S.R.T.T.U. <1396>)');
nal=animatedline('Color','r','LineWidth',1.5,'MaximumNumPoints',n);
               %ylim([-2 2]);
xlim([0 max(wtd)]);
grid on;
grid minor;
ax = gca;
ax.YLabel.String  = 'Amplitude of Bnet';
ax.YLabel.FontSize = 10;
ax.YLabel.FontWeight ='bold';
ax.XTickMode = 'manual';
ax.XTick = linspace(0,max(wtd),13);
ax.XLabel.String = 'Position of Bnet(0-360 Degree)';
ax.XLabel.FontSize = 10;
ax.XLabel.FontWeight='bold';

subplot(3,4,[1,2,5,6]);
polarplot(PBn,'r','LineWidth',2);
title('Magnetic Fields');
ax = gca;
ax.ThetaAxisUnits = 'degrees';
ax.ThetaTick = linspace(0,360,13);
ax.RLim = [0 2];
grid on;
grid minor;
hold on;

polarplot(c1 ,'g','LineWidth',1);
polarplot(c2 ,'g','LineWidth',1);
polarplot(PBn,'r','LineWidth',1);

%%

for ii=1:n        
    a(ii)=polarplot([0 theta(ii)], [0 Bn(ii)],'r','Marker','o','LineWidth',1.5,'Visible','off');    
    b(ii)=polarplot([0        0] , [0 Br(ii)],'k','Marker','o','LineWidth',1.0,'Visible','off');
    c(ii)=polarplot([0  -2*pi/3] , [0 Bs(ii)],'b','Marker','o','LineWidth',1.0,'Visible','off');
    d(ii)=polarplot([0   2*pi/3] , [0 Bt(ii)],'m','Marker','o','LineWidth',1.0,'Visible','off');
end

%%
toc;tic;
%%

for ii=1:n
    
    a(ii).Visible='on';
    b(ii).Visible='on';
    c(ii).Visible='on';
    d(ii).Visible='on';
    
    subplot(3,4,[3,4]);
    addpoints(ral,t(ii),Br(ii));
    %drawnow limitrate ;
 
    subplot(3,4,[7,8]);
    addpoints(sal,t(ii),Bs(ii));
    %drawnow limitrate ;    
    
    subplot(3,4,[11,12]);
    addpoints(tal,t(ii),Bt(ii));
   % drawnow limitrate ;

    subplot(3,4,[9,10]);
    addpoints(nal,wtd(ii),Bn(ii));
  %  drawnow limitrate ;

    if (gifc == 1)
        
        frame(ii) = getframe(1);
        im{ii} = frame2im(frame(ii));
    else
        pause(0.04);
    end
    

    
    
    %drawnow limitrate;
    
    a(ii).Visible='off';
    b(ii).Visible='off';
    c(ii).Visible='off';
    d(ii).Visible='off';
    
end
%% Movie
%      close;
%      fig=figure('NumberTitle','off','Name',['Magnetic Fields in the Stator',date],'Units','centimeters','Position',[0.25 2 40 18]);
%     movie(fig,frame,2);
%% These following codes are written for exporting plots to a .gif (animation) file

if (gifc == 1)
    %close;
%     figure;
%     for idx = 1:n
%         subplot(m,m,idx)
%         imshow(im{idx});
%     end    
    
    for idx = 1:n
        [A,map] = rgb2ind(im{idx},256);
        if idx == 1
            imwrite(A,map,filename,'gif','LoopCount',   Inf  ,'DelayTime',0);
        else
            imwrite(A,map,filename,'gif','WriteMode','append','DelayTime',0);
        end
    end
    
end
%%
close all;
toc;

%% The End... ...!?
