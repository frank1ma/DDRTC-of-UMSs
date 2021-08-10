function [out,uin,t] = rdid_nonlnr_mdl_4d(fs,tval,C,refin)
%  Date-Driven Flatness Output Searching and ADRC Control of 
%  Underactuated Nonlinear System 
%
%  Date : 07 - 12 - 2019
%  Shangjie Ma
%  --------------------------------------------------------------------
%  -Nonlinear System Case Model - 4D Minimum-phase case
%   - Parameter Def.
%   - Open-loop Matrix Def. 
%   - Stability & Controller
%   - ODE Solution

% Parameter Def.
%========================================================================
dt = 1/fs;
t = (0:dt:tval-dt)'; % time vector
N = size(t,1);

% Open-loop Matrix Def.
%========================================================================
A = [0   1  0  0;
     1   0  1  0;
     1   3  1  1;
     6   2  3  3];
B = [0;-1;0;1.5];
D = 0;
rank(ctrb(A,B))
% Stability & Controller
%========================================================================
LM_LQR_Q =diag([1 1 1 1]);
LM_LQR_R = 1;
LM_LQR_K =lqr(A,B,LM_LQR_Q,LM_LQR_R);
Ke = LM_LQR_K;

% ODE Solution
%========================================================================
x1(1) = 0;
x2(1) = 0;
x3(1) = 0;
x4(1) = 0;
u(1)= - Ke *[x1(1);x2(1);x3(1);x4(1)] + refin(1);

for i = 1:N-1
    x = [x1(i);x2(i);x3(i);x4(i)];
    x1(i+1) = x1(i) + (A(1,:)*x + B(1) * u(i) + x1(i).^3)*dt;
    x2(i+1) = x2(i) + (A(2,:)*x + B(2) * u(i))*dt;
    x3(i+1) = x3(i) + (A(3,:)*x + B(3) * u(i)+x3(i).^3)*dt;
    x4(i+1) = x4(i) + (A(4,:)*x + B(4) * u(i))*dt;
    u(i+1) = -Ke * x + refin(i+1);
end
Acl = A-B*Ke;
eigval_Acl = eig(Acl)
disp('Acl,B,C,D');
disp(Acl);
disp(B);
disp(C);
disp(D);
plot(t,x1,t,x2,t,x3,t,x4)
title('Response of System')
xlabel('time(/s)')
ylabel('response')
legend('x1','x2','x3','x4')
x = [x1;x2;x3;x4];
out = C * x;
uin = refin;




    