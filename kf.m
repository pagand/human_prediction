function [x_hat,dx_hat,fe_hat,xest,Pout]=kf(xmeas,u,xminus,Pminus)
% initilization
stdmeas=(0.005); % standard deviation of the measurement errors
stdmodel1=0.01; % standard deviation of the evolution model
stdmodel2=0.01; % standard deviation of the evolution model
dt=0.001;
M=5; C=2;K=4;
Ke=2;Be=3; ne=1.7;
a1=-C/M; a2=-K/M; b=1/M;
xest=xminus;
k11=-dt/M*ne*sign(xest(1))*(abs(xest(1)))^(ne-1)*(Ke+Be*xest(2));
k22=-dt/M*Be*sign(xest(1))*(abs(xest(1)))^ne;
F=[1 dt;a2/a1*(1-exp(a1*dt))+k11 exp(a1*dt)+k22];
H=[1 0];
B=[0; -b/a1*(1-exp(a1*dt))];
Q=[stdmodel1^2 0 ;...
0 stdmodel2^2];
R=(stdmeas)^2;
%%%%%%%%%%%%%%%%%%%%%%%
% Prediction Kalman %
%%%%%%%%%%%%%%%%%%%%%%%
xminus=F*xminus+B*u;
Pminus=F*Pminus*F'+Q;
%%%%%%%%%%%%%%%%%%%%%%
% Update Kalman %
%%%%%%%%%%%%%%%%%%%%%%
Gainm=Pminus*H'*inv(H*Pminus*H'+R); %#ok<*MINV>
xnew=xminus+Gainm*(xmeas-H*xminus);
Pnew=(eye(2)-Gainm*H)*Pminus;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Storing the predicted variables %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xest=xnew;
Pout=Pnew;
x_hat=xest(1);
dx_hat=xest(2);
fe_hat=Ke*sign(xest(1))*(abs(xest(1)))^ne+Be*sign(xest(1))*(abs(xest(1)))^ne*xest(2);
end