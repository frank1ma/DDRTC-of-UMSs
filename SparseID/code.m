clear
addpath('Functions')
addpath('Data') 
% load('alpha_exp_1k_60s_p1_use.mat')
% load('theta_exp_1k_60s_p1_use.mat')
 n=1000;
% u=theta(2,1:n);
% alpha=alpha(2,1:n);
% t=theta(1,1:n);
% theta=theta(3,1:n);
load('sin_alpha.mat')
t=alpha(1:n,1)';
alpha=alpha(1:n,2)';
load('sin_theta.mat')
theta=theta(1:n,2)';
load('sin_u.mat')
u=u(1:n,2);




% P1=cumtrapz(t,cumtrapz(t,alpha));
% P2=cumtrapz(t,cumtrapz(t,theta));
% P3=cumtrapz(t,alpha);
% P4=cumtrapz(t,theta);
F=(theta+1.039*alpha)';


P11=cumtrapz(t,cumtrapz(t,cumtrapz(t,cumtrapz(t,F))));
P12=cumtrapz(t,cumtrapz(t,cumtrapz(t,F)));
P13=cumtrapz(t,cumtrapz(t,F));
P14=cumtrapz(t,F);
Pt=cumtrapz(t,cumtrapz(t,cumtrapz(t,cumtrapz(t,u))));

% P=[P1' P2' P3' P4'];
% coeff = inv(P'*P)*P'*F
% coeff2= pinv(P)*F;
% coeff_m = -45.4377 
% 
P=[P11 P12 P13 P14 Pt];
coeff = inv(P'*P)*P'*F
%coeff2= pinv(P)*F;
%coeff_m = -45.4377 

