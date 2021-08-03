function y = fal(e,alpha,sigma)
% nonlinear feedback funtion
if abs(e) > sigma
   y = abs(e)^alpha*sign(e);
else
   y = e/(sigma^alpha);
end
