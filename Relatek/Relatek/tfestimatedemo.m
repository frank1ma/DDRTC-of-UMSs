% A= [0 1 0 0 0;
%     0 0 1 0 0;
%     0 0 0 1 0;
%     0 0 0 0 1;
%     -1 -2 -3 -4 -1];
% B = [0 ;1; 0; 1;0];
% C = [1 0 2 0 0];
% D = 0;
%[num,den] = ss2tf(A,B,C,D);
%G = tf(num,den)

% =====================================================================
% s=tf('s')
% G=((s^2+3000)*(s^2+1000))/((s^2+5*s+40)*(s^2+5*s+40)*(s^2+0.03*s+0.06))
% figure(1);
% bode(G)
% t= 0:0.001:100;
% u =  0.01*randn(1,size(t,2));
% y = lsim(G,u,t);

% =====================================================================
% y1 = x1';
% y2 = x3';
% t = alpha(1,1:60000)';
% y3 = alpha(2,1:60000)';
% uf = theta(2,1:60000)'; 
 wk = kaiser(size(t,1)/2+1,35);
% % [G1,f] = cpsd(u,y1,wk,[],size(t,2)/2+1,1000);
%  [G3,f] = cpsd(uf,y3,wk,[],size(t,1)/2+1,1000);
% % [G4,f] = cpsd(u,y3,wk,[],size(t,2)/2+1,1000);
%  [G2,~] = cpsd(uf,uf,wk,[],size(t,1)/2+1,1000);  %kaiser(size(t,2),35)
 %[Gt,f] = cpsd(refin,refin,wk,[],size(t,1)/2+1,1000); 
 [Gt,f] = pwelch(refin,wk,[],size(t,1)/2+1,1000);
 plot(log10(f),log10(Gt));
%[txy,f]=tfestimate(u,x3,wk,[],size(t,2)/2+1,1000);
%[txy,f]=tfestimate(u,y,wk,[],size(t,2),1000);
% G1 = (log10(abs((txy))));
% G11 = G1(100:14000);
% f1 = log10(f);
% f11 = f1(100:14000);
% p1 = polyfit([f11(1) f11(end)],[G11(1) G11(end)],1);
% plot(f11,G11)
% %plot(f11(1),G11(1),'rx',f11(end),G11(end),'rx')
% %line([f11(1) f11(end)],[f11(1)*p1(1)+p1(2) f11(end)*p1(1)+p1(2)])
% yy =polyval(p1,f11);
% %plot(f11,yy);
% plot(f11,G11,f11,G11-yy)
% %findpeaks(abs(G11-yy),f11,'MinPeakProminence',0.5)
% % semilogx(2*pi*f,20*log10(abs((G3./G2))),'r')
% hold on
% grid on


% [txy,w] = tfestimate(u,y,wk,[],size(t,2),1000);
% plot(log10(w),20*log10(txy));
% hold on
% faxis = log10(f(1:end));
% G01 = G1 ./G2;
% G02 = G3 ./G2;
% plot(faxis,20*log(abs(G01-0.89711*G02)));
% 
% G13 = G1 - 0.89711*G3;
% G13u = 20*log(abs(G13./G2));
% plot(log10(f),G13u);
% 
%legend('true','me')
% G1p3 = G4;
% G1p3u = log10(abs(G4./G2));
% plot(faxis,G1p3u);

%[p,S]=polyfit(log10(f(2:end)),log10(abs(G4(2:end)./G2(2:end))),1)

%plot(faxis,yy);
%line([-2 3],[-2*p(1)+p(2) 3*p(1)+p(2)])
% title('estimation of transfer function')
% hold on
% grid on
% plot(faxis,G12);
%  [seg]=divpeak(log10(abs(G1(2:20000)./G2(2:20000))),log10(f(2:20000)));
% 
% %
% plot(seg(2).faxis,seg(2).estf,'rx')
% [p1,S]=polyfit(seg(2).faxis,seg(2).estf,1)
% line([-2 3],[-2*p1(1)+p1(2) 3*p1(1)+p1(2)])
%plot(seg(1).faxis,seg(1).estf,'bo',seg(2).faxis,seg(2).estf,'rx')
%plot(seg(1).faxis,seg(1).estf,'bo',seg(2).faxis,seg(2).estf,'rx',seg(3).faxis,seg(3).estf,'gx',seg(4).faxis,seg(4).estf,'yx')