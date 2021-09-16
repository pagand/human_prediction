clear;clc;close all

load('V.mat');
phi = ttrValue_obs;
Grid = xs_whole;

xPrecision = [0.2, 0.2, 2*pi/15, 0.2];

% u
uPrecision = [0.5, 0.5];
Ua = -1:uPrecision(1):1;
Uw = -1.5:uPrecision(2):1.5;

%% DATA
load('synthetic_trajectory.mat')
w = U_synth(:,1);
a = U_synth(:,2);
X0 =  X_synth(1,:);
final_time = 20;


% load('realdata.mat')
% w = wr;
% a = ar;
% X0 =  [xr(1),yr(1),thetar(1),vr(1)];
% final_time = size(ar,2);

% Dynamic
X(1,:)  = dynamic (X0,a(1),w(1),dt);
for k=1:final_time-1
X (k+1,:) = dynamic (X(k,:),a(k+1),w(k+1),dt);
end

% % check real data
% subplot(221);plot(X(:,1)); hold on; plot(xr,'r');
% subplot(222);plot(X(:,2)); hold on; plot(yr,'r');
% subplot(223);plot(X(:,3)); hold on; plot(thetar,'r');
% subplot(224);plot(X(:,4)); hold on; plot(vr,'r');


% ignore xshow
Xshow = X;
% figure;
% scatter(X(:,1),X(:,2),'ro');
%% Goal (from PF)
% [x,y,alpha,v]
% G1 = [2,0.8,0,0];
% G2 = [0,2,0,0];
% G3 = [1,0,0,0];



horizon = 1;
N_part_g= 20; N_part_x=200;

xmin=floor(min(min(X(:,1)),min(X(:,2))));
xmax=ceil(max(max(X(:,1)),max(X(:,2))));

g_particles=xmin+rand(N_part_g,4)*(xmax-xmin);
g_particles(:,3) = rand(N_part_g,1)*2*pi-pi;
g_particles(1,:) = G_synth;
g_particles(:,4) = 0; % set goal-state velocities zero
%g_particles(1,:) = [-.07 -1.8 -3*pi/4 0];

%% heatmap meshgrid
heatmap_percision = 0.0025;
X_ = xmin:heatmap_percision:xmax;
Y_ = xmin:heatmap_percision:xmax;
[X_heatmap, Y_heatmap] = meshgrid(X_,Y_);
V_heatmap = zeros(size(X_heatmap));
%% start
len_u_comb=length(Uw)*length(Ua);
Beta = [0.1; 10];
%Beta = [10;80];
%Beta = [0.01; 100];
Gamma = [0.09 ;0.99];
%Gamma = [0.8;0.99]
%Gamma = 1;

P_Beta = (1/size(Beta,1))*ones(size(Beta,1),1)';
P_Gamma = (1/size(Gamma,1))*ones(size(Gamma,1),1)';

%% goal particle
% weight of particles for goal position
P_Goal = (1/size(g_particles,1))*ones(size(g_particles,1),1,1);

Pu=zeros(len_u_comb,final_time,size(Gamma,1),size(Beta,1),size(g_particles,1));
PXt1 = zeros(final_time,len_u_comb);
nstep_pred_states = zeros(len_u_comb,4,final_time);
nstep_pred_states_mean = zeros(final_time, 4);

g_parts_time_4plot = zeros(final_time, 4);
g_parts_satter_4plot = zeros(N_part_g,4,final_time);

XP_pf = zeros(final_time, 4);
XP_pf = X0;
xp_old_pf = (X0.*ones(N_part_x,4)+0.01*randn(N_part_x,4))';

% showing results
PBt(1,:) = P_Beta;
PGt(1,:) = P_Gamma;
PWt(1,:) = P_Goal;

% First step set Goal equal to Gg
G = g_particles;
PXp = ones(size(Pu));
% Min/Max/dx from V.mat



for k=1:final_time
    
    % prediction of u before observation
    for m=1:length(Beta)
        beta = Beta(m);
        for n=1:length(Gamma)
            gamma = Gamma(n);
            for i=1:size(G,1)
               % [Q,U] = evaluateQ(X(k,:),Uw,Ua,dt,G1,gamma);
                [Q,U] = evaluateQ2(phi,Grid,X(k,:),Uw,Ua,dt,G(i,:),gamma,horizon);
                % i > m > n :MSB
                %^P(:,k,4*(n-1)+2*(m-1)+i) = evaluateP(Q,beta); %P(:,k,n,m,i)
                Pu(:,k,n,m,i) = evaluateP(Q,beta);
            end
        end
    end
%     % for values with very small Q we set Pu=0;
%     Pu(isnan(Pu)) = 0;
%     Pu=Pu./sum(Pu,1);
    % update parameters after observing u
   % u = [w(k),a(k)]; %delete it
%     ia = find(U(:,1)==w(k));
%     ib = find(U(:,2)==a(k));
    % if the actual u is not in the control samples find nearest
    ia = find(min(abs(U(:,1)-w(k)))==abs(U(:,1)-w(k)));
    ib = find(min(abs(U(:,2)-a(k)))==abs(U(:,2)-a(k)));
    
    for t=1:length(ia)
        aa=find(ib==ia(t));
        if(~isempty(aa))
            break;
        end
    end
    ind = ib(aa);
    idx_real_u(k) = ind; 
    
    

 
    PG = zeros(length(Gamma),1);
    for gamma=1:length(Gamma)
        PG(gamma) = sum(sum(Pu(ind,k,gamma,:,:))); %+Pu(ind,k,gamma,2,1)+Pu(ind,k,gamma,1,2)+Pu(ind,k,gamma,2,2);
        % filter role
        P_Gamma(gamma) = P_Gamma(gamma)*PG(gamma);
    end
    
    % normalization
    SUM = sum(P_Gamma);
    P_Gamma = P_Gamma/SUM;
    
    for gamma=1:length(Gamma)
    
        PGt(k+1,gamma) = P_Gamma(gamma);
    end
 
    % update beta
    % Find the Q for the observed u
%     [Q1,U1] = evaluateQ2(phi,Min,Precision, X(k,:),u(1),u(2),dt,G1,Gamma(1));
%     [Q2,U2] = evaluateQ2(phi,Min,Precision, X(k,:),u(1),u(2),dt,G1,Gamma(2));
%     PP = [PGamma(1)*WG(1) PGamma(1)*WG(2) PGamma(2)*WG(1) PGamma(2)*WG(2)];
%     Pb1g1 = sum([P(ind,k,1) P(ind,k,5)].*[PGamma(1) PGamma(2)]*WG(1));
%     Pb1g2 = sum([P(ind,k,2) P(ind,k,6)].*[PGamma(1) PGamma(2)]*WG(2));
%     Pb2g1 = sum([P(ind,k,3) P(ind,k,7)].*[PGamma(1) PGamma(2)]*WG(1));
%     Pb2g2 = sum([P(ind,k,4) P(ind,k,8)].*[PGamma(1) PGamma(2)]*WG(2));
%      sm = sum(PBeta.*[Pb1g1 P1g2 Pb2g1 Pb2g2]);
    
    Pbg = zeros(length(Beta) , size(G,1));    
    for goal=1:size(G,1)
        for beta=1:length(Beta)
            for gamma=1:length(Gamma)
                Pbg(beta, goal) = Pbg(beta, goal) + Pu(ind,k,gamma,beta,goal)*P_Gamma(gamma);
            end
        end
    end
%     Pb1g1 = P(ind,k,1,1,1)+P(ind,k,2,1,1);
%     Pb1g2 = P(ind,k,1,1,2)+P(ind,k,2,1,2);
%     Pb2g1 = P(ind,k,1,2,1)+P(ind,k,2,2,1);
%     Pb2g2 = P(ind,k,1,2,2)+P(ind,k,2,2,2);
    
    for beta=1:length(Beta)
        % filtering role
        P_Beta(beta) = P_Beta(beta)*sum(Pbg(beta,:));
    end
    
%     PBeta(1) = PBeta(1)*Pb1g1;
%     PBeta(2) = PBeta(2)*Pb1g2;
%     PBeta(3) = PBeta(3)*Pb2g1;
%     PBeta(4) = PBeta(4)*Pb2g2;
    
    % normalization
    sm = sum(P_Beta);
    P_Beta = P_Beta/sm;
  
    
    % save purpose only
    for beta=1:length(Beta)
        PBt(k+1,beta) = P_Beta(beta);
    end
    
    % update g
    for goal=1:size(P_Goal,1)
        P_Goal(goal) = P_Goal(goal)*sum(Pbg(:, goal));
    end
    
    % normalize 
    sm = sum(P_Goal);
    P_Goal = P_Goal./sm;
    
%     G(:,1:2) = G(:,1:2) .* PGoal 
%     G(:,1:2) = xmin+G(:,1:2)*(xmax-xmin)
    
%     v1 = std(G(:,1),PGoal)
%     m1 = sum(PGoal.*G(:,1))/sum(PGoal)
%     pd1 = makedist('normal','mu',m1,'sigma',v1);
%     
%     v2 = std(G(:,1),PGoal*100)
%     m2 = sum((PGoal*100).*G(:,2))/sum(PGoal*100)    
%     pd2 = makedist('normal','mu',m2,'sigma',v2);
%     
%     x_values = 1:9;
%     G(:,1) = pdf(pd1,x_values)
%     G(:,2) = pdf(pd2,x_values)

    
    % save purpose only
    for goal=1:size(P_Goal,1)
        PWt(k+1,goal) = P_Goal(goal);
    end
    % new goal sampling from goal distribution
    
    % expected position in 2 step ahead (dynamic)
    % not biased :: learning rules of movement rather than dynamic
    nstep_pred_states(:,:,k) = nstepdynamic(horizon,X(k,:),U(:,2),U(:,1),dt); %TODO
    Pt = ones(size(Gamma,1),size(Beta,1),size(G,1));
    for gamma=1:size(Gamma,1)
        Pt(gamma,:,:) = Pt(gamma,:,:)*PGt(k+1,gamma);
    end    
    for beta=1:size(Beta,1)
        Pt(:,beta,:) = Pt(:,beta,:)*PBt(k+1,beta);
    end
    for goal=1:size(G,1)
        Pt(:,:,goal) = Pt(:,:,goal)*PWt(k+1,goal);
    end
    
%     PXp = ones(size(Pu));
    for jj =1:size(Pu,1)
        PXp(jj,k,:,:,:) = Pt; 
    end
    PXp(:,k,:,:,:) = Pu(:,k,:,:,:).*PXp(:,k,:,:,:);
    PXt1(k,:) =sum(sum(sum(PXp(:,k,:,:,:),3),4),5);
%     % Point estimate: average over all possibility
    nstep_pred_states_mean(k,:) = PXt1(k,:)*nstep_pred_states(:,:,k);
%% heatmap instead of point estimate taher
    points = nstep_pred_states(:,1:2,k);
    w = PXt1(k,:)';
    F = scatteredInterpolant(points,w);
    
    x = min(points(:,1)):heatmap_percision:max(points(:,1));
    y = min(points(:,2)):heatmap_percision:max(points(:,2));
    
    [Xq,Yq] = meshgrid(x,y);
%     [Xq,Yq] = meshgrid(-2:0.125:2);
    Vq = F(Xq,Yq);
    
    x0 = floor((x(1)-xmin)/heatmap_percision);
    ind_x = x0:x0+length(x)-1;
    y0 = floor((y(1)-xmin)/heatmap_percision);
    ind_y = y0:y0+length(y)-1;
    V_heatmap(ind_y, ind_x) = Vq;
    figure(1);
    surf(X_heatmap,Y_heatmap, V_heatmap);
    shading interp
    xlabel('X','fontweight','b'), ylabel('Y','fontweight','b');
    zlabel('Value - V','fontweight','b');
    view(2);
    hold on; scatter3(Xshow(1,1),Xshow(1,2),0,50,'ro','LineWidth',2);
    %%
    
    % second approach particle filter for 2step prediction
    [xestsir, stdsir, xpartires, xpartires_1step]=pf_x(X(k,:),xp_old_pf,U,PXt1(k,:),N_part_x,horizon,dt);
    xp_old_pf = xpartires_1step;
    
    XP_pf(k,:) = xestsir;
    
    % combine Gg and past goal
    
    Pg = reshape(sum(sum(PXp(ind,k,:,:,:),3),4),[N_part_g,1]);

    [xestsir,stdsir,xpartires]=pf_goal(G,Pg,g_particles,N_part_g,dt);
    %xpartires_total = xpartires;
    % choose the ng particles randomly
    
    G = xpartires;
    G(1,:) = G_synth;
    g_parts_time_4plot(k,:) = xestsir;
    g_parts_satter_4plot(:,:,k) = G;
    
    % prediction
    if mod(k,5)==1
        i = k;
        figure;
        title(i);
        hold on;plot(nstep_pred_states_mean(1:i,1),nstep_pred_states_mean(1:i,2),'b:','LineWidth',2);
        hold on;scatter(nstep_pred_states_mean(1:i,1),nstep_pred_states_mean(1:i,2),50,'bo','LineWidth',2);
        %x particles
        hold on; scatter(xpartires_1step(1,:),xpartires_1step(2,:),'gd','filled', 'LineWidth' , 10);

        [delta_x_g, delta_y_g] = pol2cart(G(:,3),0.5);
        quiver(G(:,1),G(:,2),delta_x_g,delta_y_g,0,'linewidth',2, 'MaxHeadSize',1.5)
        hold on; scatter(G(:,1),G(:,2));

        hold on;plot(Xshow(1:i,1),Xshow(1:i,2),'r','LineWidth',2);
        hold on;scatter(Xshow(1:i,1),Xshow(1:i,2),50,'ro','LineWidth',2);
        hold on;
        [delta_x, delta_y] = pol2cart(nstep_pred_states_mean(1:i,3),nstep_pred_states_mean(1:i,4)/5);
        quiver(nstep_pred_states_mean(1:i,1),nstep_pred_states_mean(1:i,2),delta_x,delta_y,0,'linewidth',2, 'MaxHeadSize',1.5)
        hold on;
        grid;
    end
%     F(k) = getframe(gcf) ;
%     drawnow
    
end

%% video saving

% % create the video writer with 1 fps
%   writerObj = VideoWriter('myVideo.avi');
%   writerObj.FrameRate = 10;
%   % set the seconds per image
% % open the video writer
% open(writerObj);
% % write the frames to the video
% for i=1:length(F)
%     % convert the image to a frame
%     frame = F(i) ;    
%     writeVideo(writerObj, frame);
% end
% % close the writer object
% close(writerObj);

%% results

time = 1:final_time;
% lables = strings(size(U,1), 1);
% for t=1:size(U,1)
%    lables(t,1) =  sprintf('%0.2f,%0.2f',U(t,1),U(t,2));
% end
% 
% p_u_a_t = PXt1';
% for t=time
%     figure;
%     plot(p_u_a_t(:,t))
%     hold on; bar(idx_real_u(t),max(p_u_a_t(:,t)),'red');
%     labels=1:25;
%     xticks(1:25)
%     xticklabels(lables)
% end

% [~,i]=max(PXt1');
% ff=U(i,:);
% subplot(211);plot(ff(:,1));
% hold on; plot(w,'r');
% legend('max estimate','real');
% ylabel('action w');
% subplot(212);plot(ff(:,2));
% hold on; plot(a,'r');
% ylabel('action a');
% xlabel('Time(Sec)');

figure;
subplot(221)
plot(g_parts_time_4plot(:,1))
ylabel('goal x')
subplot(222)
plot(g_parts_time_4plot(:,2))
ylabel('goal y')
subplot(223)
plot(g_parts_time_4plot(:,3))
ylabel('goal  \theta')
subplot(224)
plot(g_parts_time_4plot(:,4))
ylabel('goal v')

% 
figure;
subplot(211);plot(PBt(:,1))
title('Probability of Beta')
ylabel('Beta = 0.1')
subplot(212);plot(PBt(:,2))
ylabel('Beta = 10')

figure;
subplot(211);plot(PGt(:,1))
ylabel('gamma = 0.09')
title('Probability of Gamma')
subplot(212);plot(PGt(:,2))
ylabel('gamma = 0.99')

% figure;
% subplot(311);plot(PWt(:,1))
% ylabel('Goal 1 ')
% title('Goal Probability')
% subplot(312);plot(PWt(:,2))
% ylabel('Goal 2 ')
% subplot(313);plot(PWt(:,3))
% ylabel('Goal 3 ')
% subplot(914);plot(PWt(:,4))
% ylabel('Goal 4 ')
% subplot(915);plot(PWt(:,5))
% ylabel('Goal 5 ')

% subplot(916);plot(PWt(:,6))
% ylabel('Goal 6 ')
% subplot(917);plot(PWt(:,7))
% ylabel('Goal 7 ')
% subplot(918);plot(PWt(:,8))
% ylabel('Goal 8 ')
% subplot(919);plot(PWt(:,9))
% ylabel('Goal 9 ')


% h = figure;
% for t=1:Tf
%     hold on;
%     Xp1 = reshape(Xp(:,1,t),[length(Uw),length(Ua)]);
%     Xp2 = reshape(Xp(:,2,t),[length(Uw),length(Ua)]);
%     PXt = reshape(PX(t,:),[length(Uw),length(Ua)]);
%     surf(Xp1,Xp2,10*log10(PXt));
%   
%    newmap = jet;                    %starting map
% ncol = size(newmap,1);           %how big is it?
% zpos = 1 + floor(2/3 * ncol);    %2/3 of way through
% newmap(zpos,:) = [1 1 1];        %set that position to white
% colormap(newmap); 
% end

% hold on;scatter(G(1,1),G(1,2),50,'go','LineWidth',8);
% text(G(1,1),G(1,2), 'G1', 'Fontsize', 10);
% 
% for i=2:size(G,1)
% hold on;scatter(G(i,1),G(i,2),50,'yo','LineWidth',8);
% text(G(i,1),G(i,2), 'G'+string(i), 'Fontsize', 10);    
% end


 figure;
 subplot(221);
  plot(Xshow(:,1))
 hold on; plot([X0(1);X0(1);nstep_pred_states_mean(:,1)],'r:')
 hold on; plot(XP_pf(:,1),'k-.')
 legend('ground truth','weighted average','pf')
  subplot(222);
  plot(Xshow(:,2))
 hold on; plot([X0(2);X0(1);nstep_pred_states_mean(:,2)],'r:')
 hold on; plot( XP_pf(:,2),'k-.')
  subplot(223);
  plot(Xshow(:,3))
 hold on; plot([X0(3);nstep_pred_states_mean(:,3)],'r:')
 hold on; plot(XP_pf(:,3),'k-.')
  subplot(224);
 plot(Xshow(:,4))
 hold on; plot([X0(4);nstep_pred_states_mean(:,4)],'r:')
 hold on; plot(XP_pf(:,4),'k-.')

  time1 = repmat(time,[N_part_g,1]);
   time1 = reshape(time1,[1,N_part_g*final_time]);
 figure; % plot all goals in time
 subplot(311);
 scatter(time1,reshape(g_parts_satter_4plot(:,1,:),[1,N_part_g*final_time]),30,'cyan','LineWidth',0.5);
  hold on;scatter(time,(g_parts_time_4plot(:,1))',50,'k*','LineWidth',2);
  hold on;plot(time,G_synth(1)*ones(length(time)),'r-.')
  legend('Goal particle','Average','GT goal');
  ylabel('x');
 subplot(312);
 scatter(time1,reshape(g_parts_satter_4plot(:,2,:),[1,N_part_g*final_time]),30,'cyan','LineWidth',0.5);
  hold on;scatter(time,(g_parts_time_4plot(:,2))',50,'k*','LineWidth',2);
  hold on;plot(time,G_synth(2)*ones(length(time)),'r-.')
  ylabel('y');
 subplot(313);
  scatter(time1,reshape(g_parts_satter_4plot(:,3,:),[1,N_part_g*final_time]),30,'cyan','LineWidth',0.5);
  hold on;scatter(time,(g_parts_time_4plot(:,3))',50,'k*','LineWidth',2);
  hold on;plot(time,G_synth(3)*ones(length(time)),'r-.');
  ylabel('\theta');
%  subplot(224);
%  scatter(time1,reshape(g_parts_satter_4plot(:,4,:),[1,N_part_g*final_time]),30,'go','LineWidth',8);
%   hold on;scatter(time,(g_parts_time_4plot(:,4))',50,'k*','LineWidth',8);
%   ylabel('v');
  
norm(Xshow(:,1)-XP_pf(:,1))
norm(Xshow(:,2)-XP_pf(:,2))
norm(Xshow(:,3)-XP_pf(:,3))
norm(Xshow(:,4)-XP_pf(:,4))


norm(Xshow(:,1)-[X0(1);X0(1);nstep_pred_states_mean(1:end-2,1)])
norm(Xshow(:,2)-[X0(2);X0(2);nstep_pred_states_mean(1:end-2,2)])
norm(Xshow(:,3)-[X0(3);nstep_pred_states_mean(1:end-1,3)])
norm(Xshow(:,4)-[X0(4);nstep_pred_states_mean(1:end-1,4)])

