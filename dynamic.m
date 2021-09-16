function [Xout] = dynamic(X,a,w,dt)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
x = X(:,1);
y = X(:,2);
theta = X(:,3);
v = X(:,4);
xo = x + v.*cos(theta)*dt+0*a;
yo = y + v.*sin(theta)*dt+0*a;
thetao = theta+ w*dt;
vo = v+ a*dt ;
Xout = [xo,yo,thetao,vo];
end

