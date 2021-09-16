function [Xout] = nstepdynamic(horizon,X,a,w,dt)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
x = X(1);
y = X(2);
alpha = X(3);
v = X(4);
xo = x + v*cos(alpha)*dt+0*a;
yo = y + v*sin(alpha)*dt+0*a;
alphao = alpha+ w*dt;
vo = v+ a*dt ;
%for n TODO!
xo2 = xo + vo.*cos(alphao)*dt+0*a;
yo2 = yo + vo.*sin(alphao)*dt+0*a;
alphao2 = alphao+ w*dt;
vo2 = vo+ a*dt ;
Xout = [xo2,yo2,alphao2,vo2];
end

