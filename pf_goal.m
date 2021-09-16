function [xestsir,stdsir,xpartires]=pf_goal(xpartiold,prob,particle_set,n_part_g,dt)
%%
% Nparti: number of particles
% meas: meas(1), meas(2), meas(3) are the measured x,y,theta position
% xpartiold: is the Nparti of previous particles
% u: control signal u(1)=w , u(2)=a

%%
% initilization
xpartiold =xpartiold';
Nparti=n_part_g;
Nparti1=1/Nparti;
alpha = 0.25; % chance of choosing from fix particle set
stdmeas=(0.1); % standard deviation of the measurement errors


stdmodel1=0.01; % standard deviation of the evolution model
stdmodel2=0.01; % standard deviation of the evolution model
stdmodel3=0.01; % standard deviation of the evolution model
stdmodel4=0; % standard deviation of the evolution model


stdrw=0.001; % random walk of the estimated parameter (f)

% initial condition
xestsir=[0;0;0;0];stdsir=[1;1;1;1];
xpartinew=zeros(4,Nparti);
wparti=zeros(1,Nparti);
wpartin=zeros(1,Nparti);
xpartires=zeros(4,Nparti);
wpartires=zeros(1,Nparti);

%%%%%%%%%%%%%%%%%%%%%%%%
% Advancing Particles %
%%%%%%%%%%%%%%%%%%%%%%%%

for i=1:Nparti
%%%%%%%%%%%%%%%%%%%%
% 1 - PREDICTION %
%%%%%%%%%%%%%%%%%%%%
F=[1 0 0 0;0 1 0 0;0 0 1 0;0 0 0 1];

%  radndom walk
if rand>alpha
xpartinew(:,i)=F*(xpartiold(:,i))+randn*[stdmodel1;stdmodel2;stdmodel3;stdmodel4];
else
    index = randperm(size(particle_set,1),1);
    xpartinew(:,i)=particle_set(index,:);
end 
% else
%     xpartinew(:,i)=randn*[stdmodel;0];
%  end
  
%   fe(1,i)=Ke*sign(xpartinew(1,i))*(abs(xpartinew(1,i)))^(ne)+Be*sign(xpartinew(1,i))*(abs(xpartinew(1,i)))^(ne)*xpartinew(2,i)+randn*stdrw1;
%  if ( xpartinew(1,i)<0)
%      fe(1,i)=0;
%  end

      
 %%%%%%%%%%%%%%%%%%%%
% Weights %
%%%%%%%%%%%%%%%%%%%%
%wparti(i)=exp(-((xpartiold(1,i)-meas(1))/stdmeas)^2-((xpartiold(2,i)-meas(2))/stdmeas)^2-((xpartiold(3,i)-meas(3))/stdmeas)^2);
end
wparti = prob;

wtotal=sum(wparti);
wpartin=wparti./wtotal;
%%%%%%%%%%%%%%%%
% 2 - UPDATE %
%%%%%%%%%%%%%%%%
%xpartiw=xpartinew.*wparti;
%%%%%%%%%%%%%%%%%%%%%%
% Mean at time k+1 %
%%%%%%%%%%%%%%%%%%%%%%

%xestpf=xpartinew*wpartin';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Standard Deviation at time k+1 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3 - RESAMPLING %
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
stdsir(2)=std(xpartires(2,:)); % y
stdsir(3)=std(xpartires(3,:)); % theta
stdsir(4)=std(xpartires(4,:)); % v

%fe_hat=xestsir(3);
%%%%%%%%%%%%%%%%%%%%%
% End of resampling %
%%%%%%%%%%%%%%%%%%%%%
%xpartiold=xpartinew;
xestsir = xestsir';
xpartires = xpartires';
end