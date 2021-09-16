close all;
clear;

%% load grid and ttr values
dim = 4;
load('V.mat');
phi = ttrValue_obs;
Grid = xs_whole;
%%
dt = 0.2;
%%
gamma = 0.99;
beta = 10;
horizon = 1;
uPrecision = [0.5, 0.5];
Uw = -1.5:uPrecision(1):1.5;
Ua = -1:uPrecision(2):1;
%%
Tf =20;

G = [0.5, 0,pi,0]; 
for beta = 10:10
for i = 1:1
for theta=0
X = zeros(Tf,4);
Umax = zeros(Tf,2);

init_X = [2, 2, theta, 0];

X(1,:) = init_X;
Idx = 0;
for k=1:Tf
    [Q,U,Xtn] = evaluateQ2(phi,Grid,X(k,:),Uw,Ua,dt,G,gamma,horizon);
    [P] = evaluateP(Q,beta);
%     figure;
%     plot(P);
    if i==1
    % select optimal action/determinisitc
    [Pmax, Idx] = max(P,[],2);
    else
    % sample from action probability distribution
    Idx = datasample(ndgrid(1:size(P,2)),1,'weights',P);
    end

    Umax(k,:) = U(Idx,:);
    X(k+1,:) = dynamic (X(k,:),Umax(k,2),Umax(k,1),dt);

    if((X(k,1)-G(1))^2 + (X(k,2)-G(2))^2 < 0.05)...
             &&((abs(X(k,3)-G(3)) < 4*pi/30)||(abs(abs(X(k,3)-G(3))-2*pi) < 4*pi/30))...
             &&(abs(X(k,4)-G(4)) < 0.3)
        break;
    end

%     fig = figure(2);
%     
%     fig_axes = axes('Parent',fig);
%     hold(fig_axes, 'on');
%     hold on; scatter(X(1,1),X(1,2),'rx','LineWidth',2);
%     [delta_x_x, delta_y_x] = pol2cart(X(1,3),0.7);
%     quiver(X(1,1),X(1,2),delta_x_x,delta_y_x,0,'b','linewidth',1, 'MaxHeadSize',1);
%     for i=1:size(Xtn,1)
%         [delta_x, delta_y] = pol2cart(Xtn(i,3), Xtn(i,4));
%         hold on; quiver(X(k,1),X(k,2),delta_x, delta_y,0,'linewidth',0.5, 'MaxHeadSize',0.2)
%         grid on;
%         set(fig_axes, 'XLim',[-8 8], 'YLim', [-8 8]);
%     end
%     [delta_x, delta_y] = pol2cart(Xtn(Idx,3), Xtn(Idx,4));
%     hold on; scatter(X(k,1),X(k,2),'r','LineWidth',1);
%     hold on; quiver(X(k,1),X(k,2),delta_x,delta_y,0,'linewidth',4, 'MaxHeadSize',1.5)
%     hold on; scatter(G(1,1),G(1,2),'gx','LineWidth',2);    
%     [delta_x_g, delta_y_g] = pol2cart(G(:,3),0.7);
%     quiver(G(1,1),G(1,2),delta_x_g,delta_y_g,0,'g','linewidth',1, 'MaxHeadSize',1);
%     hold on; plot(X(1:k,1),X(1:k,2),'k','LineWidth',1);
%     saveas(gcf,sprintf('results/synthetic_trajectory/%03d.png',k));

end

% fig2 = figure;
% fig_axes2 = axes('Parent',fig2);
% [delta_x_g, delta_y_g] = pol2cart(G(:,3),0.5);
% quiver(G(:,1),G(:,2),delta_x_g,delta_y_g,0,'linewidth',1, 'MaxHeadSize',1);
% 
% velocity = X(k,4)
% if (X(k,4))<0.01
%     velocity = 0.5
% end
% [delta_x_g, delta_y_g] = pol2cart(X(k,3),velocity);
% quiver(X(k,1),X(k,2),delta_x_g,delta_y_g,0,'linewidth',1, 'MaxHeadSize',1);
% 
% hold on; scatter(G(1),G(2),'b');
% hold on; plot(X(1:k,1),X(1:k,2),'r','LineWidth',2);
% title(sprintf('Theta:%0.3f, k=%d',theta,k));
% grid on;
% set(fig_axes2, 'XLim',[-8 8], 'YLim', [-8 8]);

end
figure(beta);

if i ==1
scatter(X(1,1),X(1,2),'rx','LineWidth',3);
hold on; scatter(G(1,1),G(1,2),'gx','LineWidth',3); 
hold on; plot(X(1:k,1),X(1:k,2),'r','LineWidth',1);
else
hold on; plot(X(1:k,1),X(1:k,2),'Color',[1 1 1]-0.1*beta,'LineWidth',1);
end
end
end
% fig3 = figure(beta);
% fig_axes3 = axes('Parent',fig3);
% grid on;
% set(fig_axes3, 'XLim',[-8 8], 'YLim', [-8 8]);

G_synth = G;
U_synth = Umax;
X_synth = X(1:k,:);
save('synthetic_trajectory.mat', 'X_synth', 'U_synth', 'G_synth','dt');

