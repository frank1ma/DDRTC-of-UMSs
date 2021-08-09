%*      *       *       *       *       *       *       *       *       *
%
%  Date-Driven Flatness Output Searching and ADRC Control of 
%  Underactuated Nonlinear System 
%
%  Date : 05 - 29 - 2019
%  Shangjie Ma
%  --------------------------------------------------------------------
%  -Linear System Case Model 
%   - Parameter Def.
%   - Notation Def.
%   - Open-loop Matrix Def. 
%   - Stability & Controller
%   - ODE Solution
%
%*      *       *       *       *       *       *       *       *       *

% Parameter Def.
%========================================================================
eps = 1.2;
gama = 93.7811;
ita = 1.3376;
alpha2 = 4.5109;
beta = 35.6727;

% Notation Def. 
%========================================================================
LM_f1 =  beta*(alpha2+eps^2)/(alpha2*ita+eps^2*ita-eps^2);
LM_g1 =  eps*gama/(alpha2*ita+eps^2*ita-eps^2);
LM_f2 =  eps*beta/(alpha2*ita+eps^2*ita-eps^2);
LM_g2 =  ita*gama/(alpha2*ita+eps^2*ita-eps^2);

% Open-loop Matrix Def.
%========================================================================
A = [0 1 0 0; LM_f1 0 0 0; 0 0 0 1; LM_f2 0 0 0];
B = [0;LM_g1;0;LM_g2];
C = [1 0 1 0];
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
x1(1) = rand;%.1*rand;
x2(1) = rand;%.1*rand;
x3(1) = rand;%.1*rand;
x4(1) = rand;%.1*rand;
LM_X0=[x1(1) x2(1) x3(1) x4(1)];
u(1)=0;
for i = 1:N-1
    x1(i+1) = x1(i) + x2(i)*dt;
    x2(i+1) = x2(i) +( LM_f1 * x1(i) + LM_g1 * u(i))*dt;
    x3(i+1) = x3(i) + x4(i)*dt;
    x4(i+1) = x4(i) + ( LM_f2 * x1(i) + LM_g2 * u(i))*dt;
    u(i+1) = -Ke * [x1(i);x2(i);x3(i);x4(i)] + REFIN(i);
end
% 
% LM_g1
% LM_g2
Acl = A-B*Ke;
% % plot(t,x1,t,x2,t,x3,t,x4)
% % legend('x1','x2','x3','x4')
% rd=[];
%  for i = -1:0.0001:0
% [num,den]=ss2tf(Acl,B,[1 0 -0.89711 0],0);
% tf(num,den)
% bode( tf(num,den))
%grid on
%     k = 1;
%     
%     while abs(num(k)) < 1e-3
%         k = k+1;
%     end
%     rd = [rd k-1];
%     end
%  plot([-1:0.0001:0],rd)

    