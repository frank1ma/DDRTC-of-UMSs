clear all;

%% Filter Parameters
% SRV02 High-pass filter in PD control used to compute velocity
% Cutoff frequency (rad/s)
wcf_1 = 2 * pi * 50;
wcf_2 = 2 * pi * 10;
%

A=[0 0 1 0;
   0 0 0 1;
   0 80.6538 -0.9231 0;
   0 -120.9808 1.3846 0;];

B=[0;0;51.5346;-49.3846];

C=[1 0 0 0;0 1 0 0];

D=[0;0];
Q1=diag([100000 100 100 100]);
R = 100;
SYS=ss(A,B,C,D);
[K,S,CLP]=lqr(A,B,Q1,R)
