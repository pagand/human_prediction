clear all; close all;

%% xy data
load('./data/vicon_hat_3_hat_3_translation.csv')
xy_data = vicon_hat_3_hat_3_translation(600:1000,2:3);
txy = [vicon_hat_3_hat_3_translation(600:1000,1), xy_data(1:end,1:2)];

%% calculate the heading angle:
myDiff = diff(xy_data);
derived_theta = atan2(myDiff(:,2), myDiff(:,1));

txyheading = [xy_data(1:end-1,:), derived_theta];
%% velocity from xy and theta from velocity
t = txy(:,1); 
len = length(t)
for i=1:len-1
    dx = xy_data(i+1,1)-xy_data(i,1);
    dy = xy_data(i+1,2)-xy_data(i,2);
    dt = t(i+1,1)-t(i,1);
    dvx = dx/dt;
    dvy = dy/dt;
    
    v(i) = sqrt(dvx^2+dvy^2);
    thetav(i) = atan2(dvy,dvx);
end
%% angular velocity from thetav
for i=1:len-2
    dt = t(i+1)-t(i);
    dtheta = thetav(i+1) - thetav(i);
    w(i) = dtheta/dt;
end
%% acceleration from v
for i=1:len-2
	dt = t(i+1)-t(i);
    dv = v(i+1) - v(i);
    a(i) = dv/dt;
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
%