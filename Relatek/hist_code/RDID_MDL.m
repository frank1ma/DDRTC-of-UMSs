%*      *       *       *       *       *       *       *       *       *
%
%  LINEAR MODEL
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

function [OUT,u,Ke] = RDID_MDL(REFIN,C,t,dt)
% Parameter Def.
%========================================================================
N = length(t);
% Notation Def. 
%========================================================================
% Open-loop Matrix Def.
%========================================================================
A = [0   1  0  ;
     0   0  1  ;
    -1  -3 -1;];
B = [1;0;1]
% Stability & Controller
%========================================================================
LM_LQR_Q =diag([10 1 1]);
LM_LQR_R = .1;
LM_LQR_K =lqr(A,B,LM_LQR_Q,LM_LQR_R);
Ke = LM_LQR_K;
% ODE Solution
%========================================================================
x1(1) = rand;%.1*rand;
x2(1) = rand;%.1*rand;
x3(1) = rand;%.1*rand;
x4(1) = rand;%.1*rand;
u(1)= - Ke *[x1(1);x2(1);x3(1)] + REFIN(1);
for i = 1:N-1
    x = [x1(i);x2(i);x3(i)];
    x1(i+1) = x1(i) + (A(1,:)*x + B(1) * u(i))*dt;  %x2(i).^3
    x2(i+1) = x2(i) + (A(2,:)*x + B(2) * u(i))*dt;
    x3(i+1) = x3(i) + (A(3,:)*x + B(3) * u(i))*dt; %cos(x1(i)).*
%   x4(i+1) = x4(i) + (A(4,:)*x + B(4) * u(i))*dt;
    u(i+1) = -Ke * x + REFIN(i+1);
end
 Acl = A-B*Ke;
% [num,den]=ss2tf(Acl,B,[1 0 -1],0)
% num(1)
% den
% tf(num,den)
% bode(tf(num,den))
% grid on
x = [x1;x2;x3];
rd =[];
% [num,den]=ss2tf(Acl,B,[1 0 -],0);
% tf(num,den)
for i = -50:0.01:50
[num,den]=ss2tf(Acl,B,[1 0 i],0);
    k = 1;
    while num(k)==0
        k = k+1;
    end
    rd = [rd k - 1];
end
plot([ -50:0.01:50],rd)
% grid on
% plot(t,x1,t,x2,t,x3)
%  legend('x1','x2','x3')
%  hold on
OUT = C * x;




    