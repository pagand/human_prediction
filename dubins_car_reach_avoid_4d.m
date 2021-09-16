clear all; clc; close all;

% This ttr code is generated for DubinsCar3D with additive disturbances
% dx = v * cos(theta) + d
% dy = v * sin(theta) + d
% dtheta = omega
% dv = a
% Control: a and omega
% disturbance: d


global a_upper; global a_lower;
global omega_max;
global target_rad; global target_v;
global target_x; global target_y;
global target_phi;
global obs_rad; 
global obs_v_min;

a_lower = -1; a_upper = 1; 
omega_max = 1.5; 

% target setting
target_rad = 0.2   ;
target_x = 0;
target_y = 0;
target_v = 0;
target_phi = 0;

% obstacle setting
obs_rad = 0.15;
obs_v_min = -0.2;
obs_v_max = 2;

% boundaries
dim = 4;
Min = zeros(dim,1);
Max = zeros(dim,1);
Min(1) = -8;
Min(2) = -8;
Min(3) = -pi; 
Min(4) = -0.2; % we defined minimum v = -0.2 to be able to compute gradients
Max(1) = 8;
Max(2) = 8;
Max(3) = pi;
Max(4) = 2.2;

% dimension will be  161   161    31    25
dx = [0.1; 0.1; 2*pi/30; 0.1];

[xs_whole N] = gridGeneration(dim, Min, Max, dx);

Max(3) = Max(3) - dx(3);
[xs N] = gridGeneration(dim, Min, Max, dx);

% initialization
phi = 100*ones(N(1),N(2),N(3), N(4));
obs = 0*ones(N(1),N(2),N(3), N(4));

% Target init value = 0
% phi(((xs(:,:,:,1) - target_x).^2 + (xs(:,:,:,2) - target_y).^2) <= target_rad^2) = 0;
% change to circle 
flag = ((target_x - xs(:,:,:,:,1)).^2 + (target_y - xs(:,:,:,:,2)).^2) <= target_rad^2 &...
        abs(target_phi - xs(:,:,:,:,3)) <= target_rad &...
        abs(xs(:,:,:,:,4)) <= target_rad; % V constraint

phi(flag) = 0;

% Add Obstacle for v < 0
obs((xs(:,:,:,:,4) < Min(4) + obs_rad) ) = 1;
% Add Obstacle for v > 2
obs((xs(:,:,:,:,4) >  Max(4) - obs_rad) ) = 1;


% Save obstacle map
% obs_map = cat(3, obs, obs(:,:,1,:));
% save("/Users/anjianli/Desktop/robotics/project/ttr-compute/view/data/obs_map.mat", "obs_map")

% Save coordinate
% xs_whole = cat(3, xs, xs(:,:,1,:)); (TODO)
% save("/Users/anjianli/Desktop/robotics/project/ttr-compute/view/data/coordinate.mat", "xs_whole")


% %LF sweeping
% mex mex_dubins_car_reach_avoid_4d.cpp; disp('dubins car 4d reach avoid mexing done!');
% % disp('mexing has already been done!');
% 
% numIter = 10000;
% TOL = 0.5;
% 
% startTime = cputime;
% tic;
% mex_dubins_car_reach_avoid_4d(phi,xs,dx,a_upper,a_lower,omega_max,numIter,TOL);
% toc;
% mex_dubins_car_reach_avoid_4d
% endTime = cputime;
% fprintf('Total execution time %g seconds\n', endTime - startTime);


%LF sweeping
mex mexLFsweep.cpp; disp('mexing done!');

% numIter = 10000;
numIter = 50;
%TOL = eps;
TOL = 0.1;

startTime = cputime;
tic;
mexLFsweep(phi,xs,dx,a_upper,a_lower,omega_max,numIter,TOL,obs);
toc;
endTime = cputime;

ttrValue_obs = cat(3, phi, phi(:,:,1,:));
% contour(ttrValue_obs(:,:,1,3), [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20], "ShowText", "on")

contour(xs(:,:,1,10,1),xs(:,:,1,10,2),phi(:,:,1,10), "ShowText", "on")
xlabel('X')
ylabel('Y')
% save('V.mat', 'phi', 'xs')
save('V.mat', 'ttrValue_obs', 'xs_whole')

