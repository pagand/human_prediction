function [Qh,U,Xtn_out ] = evaluateQ2(V,Grid,X,Uw,Ua,dt,G,gamma,horizon)


Length_U=length(Uw)*length(Ua);


ua = meshgrid(Ua,Uw);
ua=reshape(ua',[Length_U,1]);
uw = meshgrid(Uw,Ua);
uw=reshape(uw,[Length_U,1]);

U=[uw,ua];



Xt1 = dynamic(X,ua,uw,dt);
% Xt1 = dynamic(Xt1,u1,u2,dt);
Xtn_out = Xt1;
% n-step
for i=2:horizon
   % HARD
   % call evaluateP to greedy choose the current best action
   % for time t+horizen-1
end

if horizon ==1
    Xtn = Xt1;
end

%% apply the translation and rotation due to goal change
Xtn(:,3) = Xtn(:,3)-G(3);
% check rotation
Xtn_rot = [cos(G(3)) sin(G(3));
      -sin(G(3)) cos(G(3))]*([Xtn(:,1)-G(1) Xtn(:,2)-G(2)])';
Xtn(:,1) = (Xtn_rot(1,:))';
Xtn(:,2) = (Xtn_rot(2,:))';


thet = Xtn(:,3) ;
thet(thet<-pi)=thet(thet<-pi)+2*pi;
thet(thet>pi)=thet(thet>pi)-2*pi;
Xtn(:,3) = thet;

%% discount V before (geometric series)
x1 = Grid(:,:,:,:,1);
x2 = Grid(:,:,:,:,2);
x3 = Grid(:,:,:,:,3);
x4 = Grid(:,:,:,:,4);

TTR = interpn(x1, x2, x3, x4, V, Xtn(:,1), Xtn(:,2), Xtn(:,3), Xtn(:,4), 'nearest', 100);
% TTR(TTR==0)=100;
r = -dt;
V = r*ones(Length_U,1) - gamma*TTR;

if gamma < 1
    Qh = (((1-gamma^horizon)/(1-gamma)).*V)'; 
else
    Qh = V'; 
end
% G = [0,0,0,0];
% Qh = - gamma*sqrt((Xt2(:,1)-G(1)).^2+(Xt2(:,2)-G(2)).^2)';

end

