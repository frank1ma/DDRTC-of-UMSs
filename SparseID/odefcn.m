function dxdt=odefcn(t,x)
dxdt=[x(2);-4*x(2)-3*x(1)+abs(sin(t))];
end