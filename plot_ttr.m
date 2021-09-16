close all;
load('V.mat');
dim = 4;
Min = zeros(dim,1);
Max = zeros(dim,1);
Min(1) = -4;
Min(2) = -4;
Min(3) = -pi; 
Min(4) = -0.2; % we defined minimum v = -0.2 to be able to compute gradients
Max(1) = 4;
Max(2) = 4;
Max(3) = pi;
Max(4) = 2.2;

% dimension will be 81x81x20x25
dx = [0.05; 0.05; 2*pi/20; 0.1];

D1 = Min(1):dx(1):Max(1);
D2 = Min(2):dx(2):Max(2);
D3 = Min(3):dx(3):Max(3);
D4 = Min(4):dx(4):Max(4);

% [xs N] = gridGeneration(dim, Min, Max, dx);
N = size(phi);
g = createGrid(Min, Max, N, 3);

v_var = [1, 2,  3,   4,  5,  6,  7,  8];
colors= ['r', 'r','m','y','g','c','b','k'];
 for v=6
    visSetIm(g,phi(:,:,:,:),colors(find(v_var==v)),v);
 end

return
%% X-Y
for theta=[7]
%     for v=[1,2,3,4,5,6,8,10,12,14,16,21,24]
      for v=[7]
        contour(xs(:,:,1,10,1),xs(:,:,1,10,2),phi(:,:,theta,v), "ShowText", "on")
        xlabel('X')
        ylabel('Y')
        title(sprintf('theta=%0.2f v=%0.2f',D3(theta),D4(v)))
        pause(1);
%         f = gcf;
        % Requires R2020a or later
%         exportgraphics(f,sprintf('plot_ttr/x_y/theta_%d_v_%d.png',theta,v))
    end
end
%% X-V
% for y=[1,10,20,30,40,60,80]
%     for theta=[1,2,3,4,5,6,10,15,20]
%         contour(reshape(xs(:,1,1,:,1),81,21),reshape(xs(:,1,1,:,4),81,21),reshape(phi(:,y,theta,:),81,21), "ShowText", "on")
%         xlabel('X')
%         ylabel('V')
%         title(sprintf('y=%0.2f theta=%0.2f',D2(y),D3(theta)))
%         f = gcf;
%         exportgraphics(f,sprintf('plot_ttr/x_v/y_%d_theta_%d.png',y,theta))
%     end
% end
% %% Y-V
% for x=[1,10,20,30,40,60,80]
%     for theta=[1,2,3,4,5,6,10,15,20]
%         contour(reshape(xs(1,:,1,:,2),81,21),reshape(xs(1,:,1,:,4),81,21),reshape(phi(x,:,theta,:),81,21), "ShowText", "on")
%         xlabel('Y')
%         ylabel('V')
%         title(sprintf('x=%0.2f theta=%0.2f',D1(x), D3(theta)))
%         f = gcf;
%         exportgraphics(f,sprintf('plot_ttr/y_v/x_%d_theta_%d.png',x, theta))
%     end
% end
% %% Y-Theta
% for x=[1,10,20,30,40,60,80]
%     for v=[1,2,3,4,5,6,10,15,21]
%         contour(reshape(xs(1,:,:,10,2),81,20),reshape(xs(1,:,:,10,3),81,20),reshape(phi(x,:,:,v),81,20), "ShowText", "on")
%         xlabel('Y')
%         ylabel('Theta')
%         title(sprintf('x=%0.2f v=%0.2f',D1(x),D4(v)))
%         f = gcf;
%         exportgraphics(f,sprintf('plot_ttr/y_theta/x_%d_v_%d.png',x,v))
%     end
% end
% %% X-Theta
% for y=[1,10,20,30,40,60,80]
%     for v=[1,2,3,4,5,6,10,15,21]
%         contour(reshape(xs(:,1,:,1,1),81,20),reshape(xs(:,1,:,1,3),81,20),reshape(phi(:,y,:,v),81,20), "ShowText", "on")
%         xlabel('X')
%         ylabel('Theta')
%         title(sprintf('y=%0.2f v=%0.2f',D2(y),D4(v)))
%         f = gcf;
%         exportgraphics(f,sprintf('plot_ttr/x_theta/y_%d_v_%d.png',y,v))
%     end
% end
% 
% %% Theta-V
% 
% for x=[1,10,20,30,40,60,80]
%     for y=[1,10,20,30,40,60,80]
%         contour(reshape(xs(1,1,:,:,3),20,21),reshape(xs(1,1,:,:,4),20,21),reshape(phi(x,y,:,:),20,21), "ShowText", "on")
%         xlabel('Theta')
%         ylabel('V')
%         title(sprintf('x=%0.2f y=%0.2f',D1(x),D2(y)))
%         f = gcf;
%         exportgraphics(f,sprintf('plot_ttr/theta_v/x_%d_y_%d.png',x,y))
%     end
% end