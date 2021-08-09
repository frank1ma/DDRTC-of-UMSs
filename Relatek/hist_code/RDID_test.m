%*      *       *       *       *       *       *       *       *       *
%
%  Relative Degree ID 
%  Date : 06 - 12 - 2019
%  Shangjie Ma
%  ----------------------------------------------------------------------
%   name rule needed
%
%*      *       *       *       *       *       *       *       *       *
clear 
clc

% Global Configuration
%=========================================================================
dim = 2;  % Order of ODE
Dof = 2;  
ORDF =dim*Dof;% Highest Order of F
fs = 1000/1;  % Sampling frequency 
dt = 1/fs; % sampling period
tspan = 30; % time span
t = [0:dt:tspan-dt]; % time vector
N = length(t); % time size;
WINSTART = ceil(N*0.2); % window size
WINEND = floor(N*0.9);


C = eye(3);
%C = [1 0;0 1;];
% Reference 
REFIN =  1*randn(1,size(t,2));
%REFIN =zeros(1,size(t,2));
%REFIN = cos(2*pi.*t./2.*t);
%REFIN = sin(10*t) * 10;
%REFIN(1:100)=1;
%REFIN=sin(sqrt(2)*t)+sin(sqrt(5)*t)+sin(sqrt(10)*t)+sin(sqrt(40)*t)+sin(sqrt(80)*t)+sin(sqrt(120)*t)+sin(sqrt(200)*t);
%REFIN = 1*chirp(t,0,t(end),500);% %0.01*rand(1,size(t,2));%1.*chirp(t,1,t(end),100); %0.01*randn(1,size(t,2)); %
% REFIN(10001:end)=REFIN1;
% REFIN(1:1000)=linspace(0,5,1000);
% External Model
%=========================================================================
% model selection
   %[OUT,u,Ke] = RDID_MDL(REFIN,C,t,dt);
   %[OUT,u] = RDID_MDL_DUFFING(REFIN,C,t,dt);
   DDFOS_LnrModel
   
% x1 = OUT(1,:);
% x2 = OUT(2,:);
% x3 = OUT(3,:);
u_ref = REFIN;
% x4 = OUT(4,:);
% A = [0   1  0  0;
%      0   0  1  0;
%      0   0  0  1;
%    -1   -2  -3  -4];
% eig(A)
% B = [0;0;0;10];
% eig(A-B*Ke)
% [num,den]=ss2tf(A-B*Ke,B,[1 0 0 0],0);
% tf(num,den);
% sys1=tf(num,den)
% x11=lsim(sys1,REFIN,t);
% x11=x11';

% plot(t,x1,t,x11);