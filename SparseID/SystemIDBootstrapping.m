
% Copyright 2020, All Rights Reserved
% Code by Ghazaale Leylaz
% For Paper, "An Optimal Model Identification Algorithm of 
% Nonlinear Dynamical Systems with the Algebraic Method"
% by Ghazaale Leylaz, Shangjie (Frank) Ma, and Jian-Qiao Sun

clear all
close all

%% Settings
 
MaxPol=1; % Set the max order of polynomial you want to fit
CrossedProducts=0; % set if you want to consider the cross-product terms 

% The sequence of parameters is:
% If crossProduct=0: theta,theta^2, ... ,dtheta, dtheta^2, .... ,alpha, alpha^2 ...., dalpha, dalpha^2 ....
% If crossProduct=1: theta, theta^2, ... ,dtheta, dtheta^2 ....,theta*dtheta,theta*dtheta^2,... alpha alpha^2 ....dalpha dalpha^2 ....,alpha*dalpha,alpha*dalpha^2,...


%% Import functions and data folders

addpath('Functions')
addpath('Data')

swEPSfigure
swFigSize

%% Data acquisition:

np1=1;
np2=2000;
load('alpha_exp');
load('theta_exp');
%load('alpha_exp_1k_60s_p1_use');%Alpha
%alpha_degree=data_alpha(:,2);
%alpha=alpha_degree*pi/180;
t=alpha(np1:np2,1);
alpha=alpha(np1:np2,2);

%load('theta_exp_1k_60s_p1_use');% Theta
%theta_degree=data_theta(:,3);
%theta=theta_degree*pi/180;

theta=theta(np1:np2,2);

load('data_vm'); %Vm
v=data_vm(:,2);


dt=t(2);
nt=length(t);

% compute forcee
F=theta+1.039.*alpha;

%%
ratio=0.3;
L=round(ratio*nt);
%% First Derivative estimation

%%% filtering approach
sigma_theta=0.05; % gain for theta signal
sigma_alpha=0.05; % gain for alpha signal

est_theta = nleso(theta',dt,sigma_theta);
thetaf = est_theta(1,:)';
dthetaf= est_theta(2,:)';

est_alpha = nleso(alpha',dt,sigma_alpha);
alphaf = est_alpha(1,:)';
dalphaf= est_alpha(2,:)';

%%% plot
figure()
subplot(2,2,1)
plot(t,theta,t,thetaf)
ylabel('theta')
subplot(2,2,3)
plot(t,dthetaf)
ylabel('dtheta')


subplot(2,2,2)
plot(t,alpha,t,alphaf)
ylabel('alpha')
subplot(2,2,4)
plot(t,dalphaf)
ylabel('dalpha')

print -depsc Derivatives.eps

%% Library generation, split data
% split data
ratio_tst=0.5; % 50 % for training, 50 % for cross validation (AIC calculation)
ntr=round(ratio_tst*nt);
nts=nt-ntr;

f=Al(0,2,F,t)-4*Al(1,1,F,t)+2*Al(2,0,F,t);
f_tr=f(1:ntr,1);
f_ts=f(ntr+1:end,1);

%% Model searching
% K Bootstrapping
Kb=30;

% generate Lambda
numlambda = 100;
lambdastart = -4;
lambdaend = 2;
Lambda = logspace(lambdastart,lambdaend, numlambda);


parameters=cell(MaxPol);
avg_prms=cell(MaxPol,1);
std_prms=cell(MaxPol,1);


 for polyorder=1:MaxPol
         
     clear pTheta pAlpha p1 p_tr p_ts ...
         EqPrms MSEcvLambda MSEcv_min
     
     pTheta=Lib(theta,dthetaf,t,polyorder,CrossedProducts);
     pAlpha=Lib(alpha,dalphaf,t,polyorder,CrossedProducts);
     
     % organize generated library to a neety shape 
     p=[pTheta(:,1:end),pAlpha];
     
     % split data for traing and cross validation for model selection AIC
     p_tr=p(1:ntr,:);
     p_ts=p(ntr+1:end,:);
     
    for k=1:Kb
        [P,idx]=datasample(p_tr,L);
        Pf=f_tr(idx,:);
        
        % Sparse Regression
        for i=1:numlambda
            est(:,i)=sparsifyDynamics(P,Pf,Lambda(i),1);
            MSE(i)=(p_ts*est(:,i)-f_ts)'*(p_ts*est(:,i)-f_ts)/nts;
        end
        
        % Model Selection
        [MSE_min,Index]=min(MSE);
        lambdaMin(k) = Lambda(Index);
        Estimated(:,k)=sparsifyDynamics(P,Pf,lambdaMin(k),1);
        
    end
    
        % mean values
        avg_prms=mean(Estimated,2)
        % Standard deviations
        std_prms=std(Estimated,0,2)

 end
 
% %% Saving Data
% save('SystemID.mat');
% 
% %% 2D plotting
% 
% figure
% [AIC_min,AIC_min_Index]=min(AIC);
% 
% plot(1:polyorder,AIC,'-O','LineWidth',1.5)
% xlabel('Polynomial Order')
% ylabel('$AIC$')
% hold on
% plot(AIC_min_Index,AIC_min,'or','LineWidth',3)
% grid on
% % set(gca, 'YScale', 'log')
% 
% print -depsc AIC.eps
% 
% %% Parameters
% 
% Estimated_Parameters=parameters{AIC_min_Index}




