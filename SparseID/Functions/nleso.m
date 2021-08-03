function [v]=nleso(x,h,sigma)
% Nonlinear extended state observer
% h should be at least equal to time step of original signal
% betas are the feedback gains which stabilize the nonlinear ESO.

dt = h;
% sigma = 0.05;

% sigma is a value for nonlinear function fal
% one should adjust between 0.01 to 0.5 based on the performance 

% gain should be adjust carefully. Don't change them too much if not
% necessary
% beta1=100;
% beta2=300;
beta1=100;
beta2=900;
beta3=1000;
beta4=10000;
[~,sz_v] = size(x);

% It is better to know, or approxiate some initial conditions of the expected
% derivatives to avoid oscillation. For example if x=sin(t), the v(2,1) 
% should be set as 1 since x'=cos(t), x'(0)=1. 
% The default will be zeros.
v=zeros(4,sz_v);


% nonlinear eso 
for i = 1:sz_v - 1
e = v(1,i)-x(i);
v(1,i+1) = v(1,i) + (v(2,i)-beta1 * e) * dt;
v(2,i+1) = v(2,i) + (v(3,i)-beta2 * fal(e,1/2,sigma)) * dt;
v(3,i+1) = v(3,i) + (v(4,i)-beta3 * fal(e,1/4,sigma)) * dt;
v(4,i+1) = v(4,i) + (-beta4 * fal(e,1/8,sigma)) * dt;
end

end
