
% Department of Mechanical Engineering
% University of California, Merced
% ME140 Vibration and Control
% Fall 2020
% Homework #5

% The code computes and plots forced response of a car

clear all; close all; clc;

%% System response for a selected car velocity v

% System and road parameters
m=2000; % mass(Kg)
c=1000; %damping (Ns/m)
k=200000; %stiffness (N/m)

y0=0.2; %m
L=10; %m

v = input ( ' Enter the speed of the car : ' ) ;% Select a car velocity

wn=sqrt(k/m);% natural frequency
zeta=c/(2*m*wn); % damping ratio
w=2*pi*v/L;
r=w/wn;

t=0:0.01:10; % duration

X=y0*sqrt(1+4*(zeta^2)*(r^2))/sqrt((1-(r^2))^2+(2*zeta*r)^2); % Response Amplitude

phi=atan((2*zeta*r)/(1-(r^2)));

beta=atan(-(wn^2)/(2*zeta*wn*w));

x=X*cos(w*t-beta-phi);

figure()
plot(t,x,'LineWidth',2)
xlabel('Time (s)')
ylabel('x (m)')
title(['The system response for car velocity ', num2str(v),' m/s '])

grid on

%% Finding the speed v when the forced response amplitude is the largest

clear v
r=0:0.001:4;

for j=1:length(r)
    X(j)=y0*sqrt((1+(2*zeta*r(j))^2)/((1-r(j)^2)^2+(2*zeta*r(j))^2));
end

figure()
plot(r,X,'LineWidth',2)
title('FOrced response amplitude for a range of frequency ratio r.')
xlabel('Frequency ratio r')
ylabel('X')
grid on

[X_max,index]=max(X);
r_max=r(index);
v_max=L*wn*r_max/(2*pi);

fprintf('For L=%d, the largest forced response amplitude is %7.4f m when the speed is %7.4f m/s.\n',L,X_max,v_max);

%% Displacement Transmissibility

clear zeta
zeta=[0.01,0.1,0.5,0.9];

figure()
hold on

zeta=[0.05,0.1,0.5,0.7,0.9];

for i=1:length(zeta)
    for j=1:length(r)
        Td(j)=sqrt((1+(2*zeta(i)*r(j))^2)/((1-r(j)^2)^2+(2*zeta(i)*r(j))^2));      
    end
    txt = ['zeta = ',num2str(zeta(i))];
    plot(r,Td,'LineWidth',2,'DisplayName',txt)
       
end

xlabel('Frequency Ratio r')
ylabel('Td')
title('Displacement transmissibility for different damping ratio.')
hold off
legend show
grid on


























