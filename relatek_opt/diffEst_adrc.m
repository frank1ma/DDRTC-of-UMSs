function [vf, vdot] = diffEst_adrc(t, v,r,dlt_r)

% ================================================
% differential estimation algorithm based on Active Disturbance Rejection
% Control Technique (ADRC)   developed by Jing-Qing Han (2008)
%  <<自抗干扰控制技术-估计补偿不确定因素的控制技术>>
%  Jing-Qing Han, Active Disturbance Rejection Control Technique-the
%  technique for estimating and compensating the uncertainties. National
%  Defense Industry Press, 2008.
%  %  % % %  % % %  % % %  % % %  % % %  % % %  % % 
%  solving the discrete control problem:
%  x1(k+1) = x1(k) + h*x2(k)
%  x2(k+1) = x2(k) + h*(u(k) + w(k))
%  %  % % %  % % %  % % %  % % %  % % %  % % %  % % 
%  Inputs:
%               - t: time series
%               - v: the signal needed to estimate its derivative (considered as tracking reference)
%               - h: simpling stepsize
%               - r: speed factor, i.e, |u(t)| < r in bang-bang control
%               - dlt_r: constraint condition, i.e. |w(t)| < dlt_r
%               - 
% Outputs:
%               - x1: tracking output of signal v ----->  filtered signal
%               - x2: derivative of x1 -------> differential estimation of  v(t)
% ================================================

h = t(2) - t(1);                   % sampling stepsize
r1 = r + dlt_r;   

c = sqrt(r1/r);                       % damping coefficient
h1 = 30*h;                             % a parameter used to improve performance, usually h1>=h
                                            % when h1 = h, time optimal feedback control
x = zeros(2,length(t)+1);
k = 1;
while k*h <= t(end)
    d = r*h1^2;
    a0 = h1*c*x(2,k);
    y = (x(1,k) - v(k)) + a0;
    a1 = sqrt(d*(d+8*abs(y)));
    a2 = a0 + sign(y)*(a1-d)/2;
    sy = fsg(y,-d,d);
    a = (a0 + y - a2)*sy + a2;
    sa = fsg(a,-d,d);
    fhan = -r*(a/d - sign(a))*sa - r*sign(a); 
    
    x(1,k+1) = x(1,k) + h*x(2,k);
    x(2,k+1) = x(2,k) + h*fhan;
    
    k = k+1;
end

vf = x(1,1:end-1);
vdot = x(2,1:end-1);

end 

function fsg_x = fsg(x,a,b)
fsg_x = (sign(x-a)-sign(x-b))/2;
end
