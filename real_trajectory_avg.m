clear all; close all;

%% xy data
load('./data/vicon_hat_3_hat_3_translation.csv')
xy_data = vicon_hat_3_hat_3_translation(600:1000,2:3);

ind = find(xy_data(:,1)>0.85);

dx = (.8069-.8166)/4;
dy = (.04016- .07315)/4;
xy_data(ind,1) = (.8166+dx:dx:.8166+4*dx)' ;
xy_data(ind,2) = (-.04016+dy:dy:-.04016+4*dy)' ;

% % smoothed xy data
% xy_data = smoothdata(vicon_hat_3_hat_3_translation(600:1000,2:3),1,'gaussian');

txy = [vicon_hat_3_hat_3_translation(600:1000,1), xy_data(1:end,1:2)];

%% calculate the heading angle:
myDiff = diff(xy_data);
derived_theta = atan2(myDiff(:,2), myDiff(:,1));

txyheading = [xy_data(1:end-1,:), derived_theta];
%% velocity from xy and theta from velocity
t = txy(:,1);
len = length(t);

v(1) = 0;
thetav(1) = 0;




for i=2:len-1
	dx_1 = xy_data(i,1)-xy_data(i-1,1);
    dy_1 = xy_data(i,2)-xy_data(i-1,2);
    dt_1 = t(i)-t(i-1);
    dvx_1 = dx_1/dt_1;
    dvy_1 = dy_1/dt_1;

    dx_2 = xy_data(i+1,1)-xy_data(i,1);
    dy_2 = xy_data(i+1,2)-xy_data(i,2);
    dt_2 = t(i+1)-t(i);
    dvx_2 = dx_2/dt_2;
    dvy_2 = dy_2/dt_2;
    
    dvx(i) = (dvx_1+dvx_2)/2;
    dvy(i) = (dvy_1+dvy_2)/2;
    
    v(i) = sqrt(dvx(i)^2+dvy(i)^2);
    thetav(i) = atan2(dvy(i),dvx(i));
end
%% angular velocity from thetav
w(1)=0;
for i=2:len-2
	dt_1 = t(i)-t(i-1);
    dtheta_1 = thetav(i) - thetav(i-1);
    w_1 = dtheta_1/dt_1;
    
    dt_2 = t(i+1)-t(i);
    dtheta_2 = thetav(i+1) - thetav(i);
    w_2 = dtheta_2/dt_2;
    
    w(i) = (w_1+w_2)/2;
end
%% acceleration from v
a(1) = 0;
for i=2:len-2
    dt_1 = t(i)-t(i-1);
    dv_1 = v(i) - v(i-1);
    a_1 = dv_1/dt_1;
    
    dt_2 = t(i+1)-t(i);
    dv_2 = v(i+1) - v(i);
    a_2 = dv_2/dt_2;
    
    a(i) = (a_1+a_2)/2;
end
%% plot
[x_prime,y_prime] = pol2cart(thetav(1:end), v);
xy_prime = [x_prime',y_prime'];

figure;
hold on;
scatter(xy_data(:,1),xy_data(:,2),'.');
quiver(xy_data(1:end-1,1),xy_data(1:end-1,2),xy_prime(:,1),xy_prime(:,2), 'MaxHeadSize',0.001);
xlabel('x');
ylabel('y');

figure;
plot(t(1:end-1),v)
xlabel('time');
ylabel('velocity');


figure;
plot(t(1:end-2),a)
xlabel('time');
ylabel('acc');


figure;
plot(t(1:end-2),w)
xlabel('time');
ylabel('w');

% resampling with constant time sample
t = t-t(1);
dt = 0.1;
tnew = 0:dt:3.8;
xr = interpn(t,xy_data(:,1),tnew);
yr = interpn(t,xy_data(:,2),tnew);
wr=interpn(t(1:end-2),w,tnew);
ar=interpn(t(1:end-2),a,tnew);
vr=interpn(t(1:end-1),v,tnew);
thetar =interpn(t(1:end-1), thetav,tnew);


figure;
subplot(221);plot(t,xy_data(:,1));hold on ;plot(tnew,xr,'r')
subplot(222);plot(t,xy_data(:,2));hold on ;plot(tnew,yr,'r')
subplot(223);plot(t(1:end-1),thetav);hold on ;plot(tnew,thetar,'r')
subplot(224);plot(t(1:end-1),v);hold on ;plot(tnew,vr,'r')

ar(ar>1)= 1;
ar(ar<-1)=-1;
wr(wr>1.5)=1.5;
wr(wr<-1.5)=-1.5;

figure;plot(ar);
figure;plot(wr);

save('realdata','ar','wr','dt','xr','yr','vr','thetar');