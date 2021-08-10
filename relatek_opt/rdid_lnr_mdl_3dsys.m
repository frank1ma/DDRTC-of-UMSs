function [out,uin,t] = rdid_lnr_mdl_3dsys(fs,tval,C,refin)
%  Date-Driven Flatness Output Searching and ADRC Control of 
%  Underactuated Nonlinear System 
%
%  Date : 07 - 12 - 2019
%  Shangjie Ma
%  --------------------------------------------------------------------
%  -Linear System Case Model - 3D Minimum-phase case
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
A = [0   1  0  ;
     0   0  1  ;
   -36 -36  -11;];
B = [1;0;-1];
D = 0;

% Stability & Controller
%========================================================================
Ke = [0 0 0];
% ODE Solution
%========================================================================
x1(1) = rand;
x2(1) = rand;
x3(1) = rand;
u(1)= - Ke *[x1(1);x2(1);x3(1)] + refin(1);

for i = 1:N-1
    x = [x1(i);x2(i);x3(i)];
    x1(i+1) = x1(i) + (A(1,:)*x + B(1) * u(i))*dt; 
    x2(i+1) = x2(i) + (A(2,:)*x + B(2) * u(i))*dt;
    x3(i+1) = x3(i) + (A(3,:)*x + B(3) * u(i))*dt; 
    u(i+1) = -Ke * x + refin(i+1);
end

Acl = A-B*Ke;
eigval_Acl = eig(Acl)
[num,den]=ss2tf(Acl,B,[1 0 0],0);
G = tf(num,den)
disp('Acl,B,C,D');
disp(Acl);
disp(B);
disp(C);
disp(D);
plot(t,x1,t,x2,t,x3)
title('Response of System')
xlabel('time(/s)')
ylabel('response')
legend('x1','x2','x3')
x = [x1;x2;x3];
out = C * x;
uin = refin;
end


    