
% Copyright 2020, All Rights Reserved
% Code by Ghazaale Leylaz
% For Paper, "An Optimal Model Identification Algorithm of 
% Nonlinear Dynamical Systems with the Algebraic Method"
% by Ghazaale Leylaz, Shangjie (Frank) Ma, and Jian-Qiao Sun

clear all;close all;clc;

%% Import functions and data folders

addpath('Functions')
addpath('Data')

%% Data acquisition:

np1=1;
np2=1500;

% 
% new data set
load('alpha_exp');
load('theta_exp');
t=alpha(np1:np2,1);
alpha=alpha(np1:np2,2);
theta=theta(np1:np2,2);



% load('alpha_exp_1k_60s_p1_use');%Alpha
% 
% t=alpha(1,np1:np2)';
% alpha=alpha(2,np1:np2)';
% 
% load('theta_exp_1k_60s_p1_use');% Theta
% u=theta(2,np1:np2)';
% theta=theta(3,np1:np2)';

% theta=out(1,:)';
% alpha=out(2,:)';
% u=uin;
% 
% load('sin_alpha.mat')
% sin_alpha=alpha(:,2);
% load('sin_theta.mat')
% sin_theta=theta(:,2);
% load('sin_u.mat')
% sin_input=u(:,2);
% sin_t=u(:,1);
% t=sin_t(1:np2);
% theta=sin_theta(1:np2);
% alpha=sin_alpha(1:np2);
% u=sin_input(1:np2);
% 
% 

% load('data_vm'); %Vm
% v=data_vm(:,2);

dt=1/1000;
nt=length(t);

% compute forcee
F= theta+1.037.*alpha;

%% Library generation

%flat output and its full derivatives up to 3
% P1=Al(3,3,F,t);   %F
% P2=-1.*Al(2,3,F,t)+3.*Al(3,2,F,t); %dF
% P3=-Al(1,3,F,t)+6.*Al(2,2,F,t)-6.*Al(3,1,F,t);  %ddF
% P4=-1.*Al(0,3,F,t)+9.*Al(1,2,F,t)-18.*Al(2,1,F,t)+6.*Al(3,0,F,t) %dddF
% 
% Pf=Al(3,3,theta,t);   % theta or alpha 
% P=[P1 P2 P3 P4];      % data matrix 


%flat output and its even derivatives up to 2
P1=Al(2,2,F,t);   %F
P3=2.*Al(2,0,F,t)-4.*Al(1,1,F,t)+Al(0,2,F,t); %ddF 

Pf=Al(2,2,alpha,t);   % theta or alpha 
P=[P1 P3];      % data matrix 

%flat output and its even derivatives up to 4
% P1=Al(4,4,F,t);   %F
% P2=Al(3,4,F,t)-4.*Al(4,3,F,t);%dF
% P3=1.*Al(2,4,F,t)-8.*Al(3,3,F,t)+12*Al(4,2,F,t); %ddF
% P4=Al(1,4,F,t)-12.*Al(2,3,F,t)+36.*Al(3,2,F,t)-24.*Al(4,1,F,t); %dddF
% P5=--Al(4,4,u,t);   % input
% 
% Pf=24.*Al(4,0,F,t)-96.*Al(3,1,F,t)+72.*Al(2,2,F,t)-16.*Al(1,3,F,t)+Al(0,4,F,t); %ddddF
% 
% 
% P=[P1  P3  P5];      % data matrix 



%%
ratio=0.5; % bootstrapping sample
L=round(ratio*nt);

%%  split data
% split data
ratio_tst=0.5; % 60 % for training, 40 % for cross validation 
ntr=round(ratio_tst*nt);
nts=nt-ntr;

Pf_tr=Pf(1:ntr,1);
Pf_ts=Pf(ntr+1:end,1);

P_tr=P(1:ntr,:);
P_ts=P(ntr+1:end,:);

%% Model searching
% K Bootstrapping
Kb=10;

% generate Lambda
numlambda = 100;
lambdastart = -20;
lambdaend = 2;
Lambda = logspace(lambdastart,lambdaend, numlambda);


for k=1:Kb
    
    [p_b,idx]=datasample(P_tr,L);
    pf_b=Pf_tr(idx,:);
        
    % Sparse Regression
    for i=1:numlambda
        est(:,i)=sparsifyDynamics(p_b,pf_b,Lambda(i),1);
        MSE(i)=(P_ts*est(:,i)-Pf_ts)'*(P_ts*est(:,i)-Pf_ts)/nts;
    end
        
    % Model Selection
    [MSE_min,Index]=min(MSE);
    lambdaMin(k) = Lambda(Index);
    Estimated(:,k)=sparsifyDynamics(p_b,pf_b,lambdaMin(k),1);
        
end
%     subplot(211)
%     semilogx(Lambda,est);
%     subplot(212)
%     loglog(Lambda,MSE)
% mean values
avg_prms=mean(Estimated,2)
% Standard deviations
std_prms=std(Estimated,0,2)






