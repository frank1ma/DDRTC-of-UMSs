function [out,uin,t] = rdid_lnr_mdl_ivp(fs,tval,C,refin)
%  Date-Driven Flatness Output Searching and ADRC Control of 
%  Underactuated Nonlinear System 
%
%  Date : 05 - 29 - 2019
%  Shangjie Ma
%  --------------------------------------------------------------------
%  -Linear System Case Model - Inverted Pendulum
%   - Parameter Def.
%   - Notation Def.
%   - Open-loop Matrix Def. 
%   - Stability & Controller
%   - ODE Solution
%
%  Reference:
%  Linear active disturbance rejection control of underactuated systems: 
%  The case of the Furuta pendulum
%  M.Ram¨ªrez-Neriaa1H.Sira-Ram¨ªreza1R.Garrido-Moctezumaa1A.Luviano-Ju¨¢rezb


% Parameter Def.
%========================================================================
eps = 1.2;
gama = 93.7811;
ita = 1.3376;
alpha2 = 4.5109;
beta = 35.6727;
dt = 1/fs;
t = (0:dt:tval-dt)'; % time vector
N = size(t,1);

% Notation Def. 
%========================================================================
LM_f1 =  beta*(alpha2+eps^2)/(alpha2*ita+eps^2*ita-eps^2);
LM_g2 =  eps*gama/(alpha2*ita+eps^2*ita-eps^2);
LM_f2 =  eps*beta/(alpha2*ita+eps^2*ita-eps^2);
LM_g1 =  ita*gama/(alpha2*ita+eps^2*ita-eps^2);

% Open-loop Matrix Def.
%========================================================================
A = [0 1 0 0; LM_f1 0 0 0; 0 0 0 1; LM_f2 0 0 0];
B = [0;LM_g1;0;LM_g2];
%C = [1 0 1 0];
D = 0;

% Stability & Controller
%========================================================================
LM_LQR_Q =diag([1 1 1 1]);
LM_LQR_R = 1;
LM_LQR_K =lqr(A,B,LM_LQR_Q,LM_LQR_R);
if any(eig(A)>0)
    Ke = LM_LQR_K;
else
    Ke = 0;
end
% ODE Solution
%========================================================================
x1(1) = rand;
x2(1) = rand;
x3(1) = rand;
x4(1) = rand;
u(1) = -Ke * [x1(1);x2(1);x3(1);x4(1)] + refin(1);
for i = 1:N-1
    x1(i+1) = x1(i) + x2(i)*dt;
    x2(i+1) = x2(i) +( LM_f1 * x1(i) + LM_g1 * u(i))*dt;
    x3(i+1) = x3(i) + x4(i)*dt;
    x4(i+1) = x4(i) + ( LM_f2 * x1(i) + LM_g2 * u(i))*dt;
    u(i+1) = -Ke * [x1(i);x2(i);x3(i);x4(i)] + refin(i);
end
Ke
Acl = A-B*Ke;
eig(Acl)
disp('Acl,B,C,D');
disp(Acl)
disp(B);
disp(C)
disp(D)
x = [x1;x2;x3;x4];
figure(1)
plot(t,x1,t,x2,t,x3,t,x4)
title('Response of System')
xlabel('time(/s)')
ylabel('response')
legend('x1','x2','x3','x4')
out = C * x;
uin = refin;
end
    