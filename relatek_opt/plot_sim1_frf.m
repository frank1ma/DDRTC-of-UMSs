clear;
load('_n095seg_del_zp','seg_zp','seg_del')
load('_n095segout.mat','seg','max_seg')

%% Coefficient = -0.95
f=figure;
% mark unused portion
plot(10.^(seg_del.x),10.^(seg_del.y),'-g')
hold on
% mark peaks
plot(10.^(seg_zp.x),10.^(seg_zp.y),'kx')
size_seg=size(seg,2);
color_list = colormap(lines(size_seg));
% plot segs
for k = 1:size_seg
    if(k==max_seg)
        plot(10.^(seg(k).faxis),10.^(seg(k).estf),'r-','LineWidth',3);
        continue;
    end
    plot(10.^(seg(k).faxis),10.^(seg(k).estf),'-.','LineWidth',3,'Color',color_list(k,:));
end

% axes
fax = f.CurrentAxes;
set(fax,'YLim',[10^-6.5,10^0.8])
set(fax,'XLim',[0.3,1000])
fax.XScale='log'
fax.YScale='log'
xticks([10^0 10^1 10^2])
yticks([10^-6  10^-4  10^-2 10^0 10^1])
xlabel('Frequency(Hz)')
ylabel('Magnitude')

% grid 
fax.XAxis.MinorTickValues = 0; 
fax.YAxis.MinorTickValues = 0;
grid minor;

%Legend
legend();
fax.Legend.String{1}='unused data'
fax.Legend.String{2}='break point'
fax.Legend.String{3}='segment 1 for est.'
fax.Legend.String{4}='segment 2 for est.'
%
str1 = ['$(' num2str(seg(1).ptg*100) '\%)$'];
str2 = ['$-' num2str(seg(2).slope_of_seg) '$' '/dec.' newline '$(' num2str(seg(2).ptg*100) '\%)$'];
annotation('textbox',[0.3 0.6 0.2 0.2],'String',str1,'FitBoxToText','on','LineStyle','none','Interpreter','latex','FontSize',14)
annotation('textbox',[0.65 0.25 0.2 0.2],'String',str2,'FitBoxToText','on','LineStyle','none','Interpreter','latex','FontSize',14)

%% Coefficient = -0.8971
clear;
load('_n08971seg_del_zp','seg_zp','seg_del')
load('_n08971segout.mat','seg','max_seg')
f=figure;
% mark unused portion
plot(10.^(seg_del.x),10.^(seg_del.y),'-g')
hold on
% mark peaks
plot(10.^(seg_zp.x),10.^(seg_zp.y),'kx')
size_seg=size(seg,2);
% plot segs
for k = 1:size_seg
    if(k==max_seg)
        plot(10.^(seg(k).faxis),10.^(seg(k).estf),'r-','LineWidth',3);
        continue;
    end
    plot(10.^(seg(k).faxis),10.^(seg(k).estf),'-.','LineWidth',3,'Color',rand(1,3));
end

% axes
fax = f.CurrentAxes;
set(fax,'YLim',[10^-11,10^0.8])
set(fax,'XLim',[0.3,1000])
fax.XScale='log'
fax.YScale='log'
xticks([10^0 10^1 10^2 10^3])
yticks([10^-10 10^-8 10^-6  10^-4  10^-2 10^0 10^1])
xlabel('Frequency(Hz)')
ylabel('Magnitude')

% grid 
fax.XAxis.MinorTickValues = 0; 
fax.YAxis.MinorTickValues = 0;
grid minor;

%Legend
legend();
fax.Legend.String{1}='unused data'
fax.Legend.String{2}='break point'
fax.Legend.String{3}='segment 1 for est.'
fax.Legend.String{4}='segment 2 for est.'
%
str1 = ['$-' num2str(seg(1).slope_of_seg) '$' '/dec.' newline '$(' num2str(seg(1).ptg*100) '\%)$'];
str2 = ['$-' num2str(seg(2).slope_of_seg) '$' '/dec.' '$(' num2str(seg(2).ptg*100) '\%)$'];
annotation('textbox',[0.2 0.45 0.2 0.2],'String',str1,'FitBoxToText','on','LineStyle','none','Interpreter','latex','FontSize',14)
annotation('textbox',[0.52 0.05 0.2 0.2],'String',str2,'FitBoxToText','on','LineStyle','none','Interpreter','latex','FontSize',14)

%% Coefficient = 0.35
clear;
load('_035seg_del_zp','seg_zp','seg_del','o_faxis','o_estf')
load('_035segout.mat','seg','max_seg')
f=figure;
% mark unused portion
plot(10.^(o_faxis),10.^(o_estf),'-g')
hold on
% mark peaks
plot(10.^(seg_zp.x),10.^(seg_zp.y),'kx')
size_seg=size(seg,2);
color_list = colormap(lines(size_seg));
marker_list = ['-.','*','--','-'];
% plot segs
for k = 1:size_seg
    if(k==max_seg)
        plot(10.^(seg(k).faxis),10.^(seg(k).estf),'r-','LineWidth',3);
        continue;
    end
    plot(10.^(seg(k).faxis),10.^(seg(k).estf),marker_list(k),'LineWidth',3,'Color',color_list(k,:));
end

% axes
fax = f.CurrentAxes;
set(fax,'YLim',[10^-6,10^1])
set(fax,'XLim',[0.3,500])
fax.XScale='log'
fax.YScale='log'
xticks([10^0 10^1 10^2 10^3])
yticks([10^-6  10^-4  10^-2 10^0 10^1])
xlabel('Frequency(Hz)')
ylabel('Magnitude')

% grid 
fax.XAxis.MinorTickValues = 0; 
fax.YAxis.MinorTickValues = 0;
grid minor;

%Legend
legend();
fax.Legend.String{1}='unused data'
fax.Legend.String{2}='segment 1.'

%
str1 = ['$-' num2str(seg(1).slope_of_seg) '$' '/dec.' '$(' num2str(seg(1).ptg*100) '\%)$'];
annotation('textbox',[0.5 0.08 0.2 0.2],'String',str1,'FitBoxToText','on','LineStyle','none','Interpreter','latex','FontSize',14)


%% Coefficient = 0.75
clear;
load('_075seg_del_zp','seg_zp','seg_del')
load('_075segout.mat','seg','max_seg')
f=figure;
% mark unused portion
plot(10.^(seg_del.x),10.^(seg_del.y),'-g')
hold on
% mark peaks
plot(10.^(seg_zp.x),10.^(seg_zp.y),'kx')
size_seg=size(seg,2);
color_list = colormap(lines(size_seg));
marker_list = ['o','*','>','-'];
% plot segs
for k = 1:size_seg
    if(k==max_seg)
        plot(10.^(seg(k).faxis),10.^(seg(k).estf),'r-','LineWidth',4);
        continue;
    end
    plot(10.^(seg(k).faxis),10.^(seg(k).estf),marker_list(k),'LineWidth',1,'Color',color_list(k,:));
end

% axes
fax = f.CurrentAxes;
set(fax,'YLim',[10^-6,10^0.8]);
set(fax,'XLim',[0.3,1000]);
fax.XScale='log';
fax.YScale='log';
xticks([10^0 10^1 10^2 10^3]);
yticks([10^-6  10^-4  10^-2 10^0 10^1]);
xlabel('Frequency(Hz)');
ylabel('Magnitude');

% grid 
fax.XAxis.MinorTickValues = 0; 
fax.YAxis.MinorTickValues = 0;
grid minor;

%Legend
legend();
fax.Legend.String{1}='unused data'
fax.Legend.String{2}='break point'
fax.Legend.String{3}='segment 1 for est.'
fax.Legend.String{4}='segment 2 for est.'
fax.Legend.String{5}='segment 3 for est.'
%
str1 = ['$(' num2str(seg(1).ptg*100) '\%)$'];
str2 = ['$(' num2str(seg(2).ptg*100) '\%)$'];
str3 = ['$-' num2str(seg(3).slope_of_seg) '$' '/dec.' '$(' num2str(seg(3).ptg*100) '\%)$'];
annotation('textbox',[0.2 0.7 0.2 0.2],'String',str1,'FitBoxToText','on','LineStyle','none','Interpreter','latex','FontSize',14)
annotation('textbox',[0.2 0.45 0.2 0.2],'String',str2,'FitBoxToText','on','LineStyle','none','Interpreter','latex','FontSize',14)
annotation('textbox',[0.53 0.08 0.2 0.2],'String',str3,'FitBoxToText','on','LineStyle','none','Interpreter','latex','FontSize',14)

swEPSfigure
swFigSize