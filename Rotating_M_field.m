close all
clear
clc
f=50;
t=0:1/(f*50):1/f;
w=2*pi*f
Bm=1;
Br=Bm*sin(w*t);
Bs=Bm*sin(w*t-2*pi/3);
Bt=Bm*sin(w*t+2*pi/3);
[Rx,Ry]=pol2cart(0,Br);
[Sx,Sy]=pol2cart(-2*pi/3,Bs);
[Tx,Ty]=pol2cart(2*pi/3,Bt);
Nx=Rx+Sx+Tx;
Ny=Ry+Sy+Ty;
[theta,Bn]=cart2pol(Nx,Ny);
c1=sin(w*t)+j*cos(w*t)
c2=1.5*c1

for i=1:length(t)
  plot(c1,'c','LineWidth',0.6);
  hold on;
  plot(c2,'g','LineWidth',0.6)
  hold on;
  %title(["Rotating M.Field By H.Najafi"]);
  %legend([pr ps pt,pn],'R','S','T','Bnet');
  plot([0 Rx(i)],[0 Ry(i)],'k','Marker','s','LineWidth',2);
  hold on;
  plot([0 Sx(i)],[0 Sy(i)],'b','Marker','s','LineWidth',2);
  hold on;
  plot([0 Tx(i)],[0 Ty(i)],'m','Marker','s','LineWidth',2);
  hold on;
  plot([0 Nx(i)],[0 Ny(i)],'r','Marker','s','LineWidth',3);
  axis square;
  axis([-2 2 -2 2]);
  %legend('r=1','r=1.5','R','S','T','Bnet');
  drawnow
  frame = getframe(1);
    im{i} = frame2im(frame);
  hold off;
end

% close;
% figure;
% nImages=length(t);
% for idx = 1:nImages
%     subplot(8,8,idx)
%     imshow(im{idx});
% end

% filename = 'testAnimated.gif'; % Specify the output file name
% for idx = 1:nImages
%     [A,map] = rgb2ind(im{idx},256);
%     if idx == 1
%         imwrite(A,map,filename,'gif','LoopCount',Inf,'DelayTime',1);
%     else
%         imwrite(A,map,filename,'gif','WriteMode','append','DelayTime',1);
%     end
% end