function [p] = Lib(x,dx,t,polyorder,CrossedProducts) 

    
for i=1:polyorder
    px(:,i)=Al(2,2,x.^i,t);
end
    
pv(:,1)=Al(1,2,x,t)-2*Al(2,1,x,t);
    
for j=2:polyorder
    pv(:,j)=Al(2,2,dx.^j,t);
end
    
for k=1:polyorder
    for l=1:polyorder
        pxv(:,(k-1)*polyorder+l)=Al(2,2,(x.^k).*(dx.^l),t);
    end
end

if CrossedProducts==1
    p=[px,pv,pxv]; 
elseif CrossedProducts==0
    p=[px,pv];
end



    