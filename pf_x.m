function [xestsir,stdsir,xpartires,xpartires_1step]=pf_x(X,xpartiold,U,prob, n_part_x, horizon,dt)
%%
% Nparti: number of particles
% meas: meas(1), meas(2) are the measured x,y position
% xpartiold: is the Nparti of previous particles
% u: control signal u(1)=w , u(2)=a

%%
% initilization

Nparti=n_part_x; 
Nparti1=1/Nparti;
std_meas=(2e-1); % standard deviation of the measurement errors

stdmodel1=0.001; % standard deviation of the evolution model
stdmodel2=0.001; % standard deviation of the evolution model
stdmodel3=0.001; % standard deviation of the evolution model
stdmodel4=0.001; % standard deviation of the evolution model



% initial condition
xestsir=[0;0;0;0];stdsir=[1;1;1;1];
xpartinew=zeros(4,Nparti);
x_onestep=zeros(4,Nparti);
wparti=zeros(1,Nparti);
wpartin=zeros(1,Nparti);
xpartires=zeros(4,Nparti);
xpartires_1step=zeros(4,Nparti);
wpartires=zeros(1,Nparti);

% sampling from actions
sum_probs = zeros(size(prob));
for i=1:size(prob,2)
    sum_probs(i) = sum(prob(1:i));
end
 % generating measurement input
 meas = X+randn(1,4)*std_meas;

%%%%%%%%%%%%%%%%%%%%%%%%
% Advancing Particles %
%%%%%%%%%%%%%%%%%%%%%%%%

for i=1:Nparti
%%%%%%%%%%%%%%%%%%%%
% 1 - ESTIMATION %
%%%%%%%%%%%%%%%%%%%%

F=[1 0 0 dt*cos(xpartiold(3,i));0 1 0 dt*sin(xpartiold(3,i)); 0 0 1 0; 0 0 0 1];

%F=[1 dt 0;0 1 0;0 0 1];
% two step prediction
rand_idx = rand;
sampeled_u_idx = find(sum_probs>rand_idx);
sampeled_u_idx = sampeled_u_idx(1);
sampeled_u = U(sampeled_u_idx,:);

x_onestep(:,i)=F*xpartiold(:,i)+[0;0;dt*sampeled_u(1);dt*sampeled_u(2)] +...
    randn*[stdmodel1;stdmodel2;stdmodel3;stdmodel4];

%%%%%%%%%%%%%%%%%%%%
% 2 - PREDICTION %
%%%%%%%%%%%%%%%%%%%%
% i =1
xpartinew(:,i) = x_onestep(:,i);
for j=2:horizon
    if j == horizon
        xpartinew(:,i) = F*xpartinew(:,i)+[0;0;dt*sampeled_u(1);dt*sampeled_u(2)]+...
                        randn*[stdmodel1;stdmodel2;stdmodel3;stdmodel4];
    else
        xpartinew(:,i) = F*xpartinew(:,i)+[0;0;dt*sampeled_u(1);dt*sampeled_u(2)];
    end
end
% else
%     xpartinew(:,i)=randn*[stdmodel;0];
%  end
  
%   fe(1,i)=Ke*sign(xpartinew(1,i))*(abs(xpartinew(1,i)))^(ne)+Be*sign(xpartinew(1,i))*(abs(xpartinew(1,i)))^(ne)*xpartinew(2,i)+randn*stdrw1;
%  if ( xpartinew(1,i)<0)
%      fe(1,i)=0;
%  end

      

% integrating the probability of each action with the stdmeas
% std_meas = std_meas/prob(i);
 %%%%%%%%%%%%%%%%%%%%
% Weights %
%%%%%%%%%%%%%%%%%%%%
wparti(i)=exp(-((x_onestep(1,i)-meas(1))/std_meas)^2-((x_onestep(2,i)-meas(2))/std_meas)^2 ...
             -((x_onestep(3,i)-meas(3))/std_meas)^2 -((x_onestep(4,i)-meas(4))/std_meas)^2);
        

end


wtotal=sum(wparti);
wpartin=wparti./wtotal;
%%%%%%%%%%%%%%%%
% 3 - UPDATE %
%%%%%%%%%%%%%%%%
%xpartiw=xpartinew.*wparti;
%%%%%%%%%%%%%%%%%%%%%%
% Mean at time k+1 %
%%%%%%%%%%%%%%%%%%%%%%

%xestpf=xpartinew*wpartin';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Standard Deviation at time k+1 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 4- RESAMPLING %
%%%%%%%%%%%%%%%%%%%%%%
%1
cresa=zeros(Nparti,1);
uresa=zeros(Nparti,1);
cresa(1)=wpartin(1);
for i=2:Nparti
cresa(i)=cresa(i-1)+wpartin(i);
end
iresa=1;
uresa(1)=rand*Nparti1;
for j=1:Nparti
uresa(j)=uresa(1)+Nparti1*(j-1);
while uresa(j) > cresa(iresa)
iresa=iresa+1;
end
xpartires(:,j)=xpartinew(:,iresa);
xpartires_1step(:,j)=x_onestep(:,iresa);
wpartires(j)=Nparti1;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Means and standard deviation %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xestsir(1)=mean(xpartires(1,:)); % x
xestsir(2)=mean(xpartires(2,:)); % y
xestsir(3)=mean(xpartires(3,:)); % theta
xestsir(4)=mean(xpartires(4,:)); % v



stdsir(1)=std(xpartires(1,:)); % x
stdsir(2)=std(xpartires(2,:)); % dx
stdsir(3)=std(xpartires(3,:)); % theta
stdsir(4)=std(xpartires(4,:)); % v


%fe_hat=xestsir(3);
%%%%%%%%%%%%%%%%%%%%%
% End of resampling %
%%%%%%%%%%%%%%%%%%%%%
%xpartiold=xpartinew;
end