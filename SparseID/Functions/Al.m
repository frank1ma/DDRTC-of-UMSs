
% this function is to inverse laplace transform
% L^-1((1/s^m)*(d^n x(s)/ds^n)=m-integral((t^n)*x)

function yout=Al(m,n,x,t)

if m==0
    yout=(t.^n).*x;
    
elseif m>=0
    
    dy(:,1)=(t.^n).*x;  
    
    for i=1:m

        %%%  Trapezoids integration
        y=cumtrapz(t,dy);
        
        % go to the next integration
        dy=y;
    end
    yout=y;

end

end


