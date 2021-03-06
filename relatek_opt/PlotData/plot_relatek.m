% Plot for Relatek
% Simulation 1 - Inverted Pendulum 
% Shangjie(Frank) Ma
% 10-24-2019


%
load('sim1_result')
load('sim1_combval')
load('sim1_radval')
load('sim1_ptgval')

% format setting
swEPSfigure
swFigSize
TitleFlag = 0;   % 1 to show title, 0 hide.
GridFlag  = 1; % 1  to turn grid on, 0 off.

% find best ratio and its max relative degree in the charts
maxrd = max(rdval);             % max rd
minrd = min(rdval);             
num_of_max_rd = find(rdval==maxrd); % possible ratio for max rd
 ptgcfdt= ptgval(num_of_max_rd)./sum(ptgval(num_of_max_rd));
 idx_of_max_ptgcfdt = find(ptgcfdt==max(ptgcfdt));
         
target_list = combval(num_of_max_rd(1) + idx_of_max_ptgcfdt - 1);
[~,target_idx] = min(abs(target_list));
target_ratio = target_list(target_idx);
target_ratio
maxrd
%% sim1 Closed-loop Response 
f1=figure(1);
% sub 1
f1Ax1=subplot(211);
fSim1ClosedloopResponseX1=plot(t,out(1,:));
hold on;
plot(t(1),out(1,1),'rx')
legend('response','initial value','Orientation','horizontal')
hold off;
ylabel('$\theta_{\delta}(rad)$')
xlabel('time$(s)$')
set(f1Ax1,'YLim',[-0.5 1]);      % axes
set(f1Ax1,'XLim',[0 20]);      % axes
if (TitleFlag)                   % title
    title('Simulation 1 Closed-loop Response-$\theta_{\delta}$');
end
if (GridFlag)                    % grid 
    grid on
end
% sub 2 
f1Ax2=subplot(212);
fSim1ClosedloopResponseX3=plot(t,out(2,:));
hold on;
plot(t(1),out(2,1),'rx')
legend('response','initial value','Orientation','horizontal')
hold off;
ylabel('$\phi_{\delta}(rad)$')
xlabel('time$(s)$')
set(f1Ax2,'YLim',[-1.5 1]);
set(f1Ax2,'XLim',[0 20]);      % axes
if (TitleFlag)
    title('Simulation 1 Closed-loop Response-$\phi_{\delta}$');
end
if (GridFlag)
    grid on
end

%% sim1 relative degree global chart
f2=figure(2);
plot(combval,rdval)
ylabel('Relative Degree($r$)')
xlabel('Ratio($C_r$)')
% axes
f2ax = f2.CurrentAxes;
set(f2ax,'YLim',[1 5]);   
% title
if (TitleFlag)            
    title('Relative Degree');
end
% grid 
if (GridFlag)             
    grid on
end 
% arrow position
 arrow_start_cor = [target_ratio+0.2 maxrd+0.3];   
 arrow_end_cor =[target_ratio+0.02 maxrd+0.05];
 width_f2ax = f2ax.Position(3);
 height_f2ax = f2ax.Position(4);
 end_point_x = ((arrow_end_cor(1)-f2ax.XLim(1))/(f2ax.XLim(2)-f2ax.XLim(1)))*width_f2ax + f2ax.Position(1);
 end_point_y = ((arrow_end_cor(2)-f2ax.YLim(1))/(f2ax.YLim(2)-f2ax.YLim(1)))*height_f2ax + f2ax.Position(2);
 start_point_x = ((arrow_start_cor(1)-f2ax.XLim(1))/(f2ax.XLim(2)-f2ax.XLim(1)))*width_f2ax + f2ax.Position(1);
 start_point_y = ((arrow_start_cor(2)-f2ax.YLim(1))/(f2ax.YLim(2)-f2ax.YLim(1)))*height_f2ax + f2ax.Position(2);
 % text
f2text = ['$r^{*}=',num2str(maxrd),'$' newline '$C_r^{*}=' num2str(target_ratio) '$'];  
% draw arrow
annotation('arrow',[ start_point_x,end_point_x],[start_point_y,end_point_y]);           
% add text
annotation('textbox',[ start_point_x,start_point_y,0.1,0.1],'String',f2text,'LineStyle','none','Interpreter','latex','FitBoxToText','on'); 

%% sim1 relative degree zoom-in chart
f3 = figure(3);
zoom_in_range = 4;
size_max_rd = size(num_of_max_rd,1);
% left
yyaxis left
stairs(combval(num_of_max_rd(1)-zoom_in_range:num_of_max_rd(size_max_rd)+zoom_in_range),rdval(num_of_max_rd(1)-zoom_in_range:num_of_max_rd(size_max_rd)+zoom_in_range),...
    '-','LineWidth',1.5)       
ylabel('Relative Degree')
xlabel('Ratio($C_r$)')
hold on
plot(target_ratio,maxrd,'ro','MarkerSize',6)
hold off
% right
yyaxis right
bar(combval(num_of_max_rd),ptgcfdt)
ylabel('Relative Data Usage')
if (TitleFlag)
    title('Relative Degree and Data Usage Zoomed In');
end
if (GridFlag)
    grid on
end
% axes
f3ax = f3.CurrentAxes;
f3ax.YAxis(1).Limits=[2 5];
f3ax.YAxis(2).Limits=[0 0.20];
f3ax.XAxis.Limits=[-0.8981 -0.8940];
%legend
legend('Relative Degree','Optimal point','Relative Data Usage','Orientation','vertical')
% arrow position
 arrow_start_cor = [target_ratio-0.0002 maxrd+0.35];   
 arrow_end_cor =[target_ratio maxrd+0.08];
 width_f3ax = f3ax.Position(3);
 height_f3ax = f3ax.Position(4);
 end_point_x = (arrow_end_cor(1)-f3ax.XLim(1))/(f3ax.XLim(2)-f3ax.XLim(1))*width_f3ax + f3ax.Position(1);
 end_point_y = ((arrow_end_cor(2)-f3ax.YAxis(1).Limits(1))/(f3ax.YAxis(1).Limits(2)-f3ax.YAxis(1).Limits(1)))*height_f3ax + f3ax.Position(2);
 start_point_x = ((arrow_start_cor(1)-f3ax.XLim(1))/(f3ax.XLim(2)-f3ax.XLim(1)))*width_f3ax + f3ax.Position(1);
 start_point_y = ((arrow_start_cor(2)-f3ax.YAxis(1).Limits(1))/(f3ax.YAxis(1).Limits(2)-f3ax.YAxis(1).Limits(1)))*height_f3ax + f3ax.Position(2);
 % text
f2text = ['$r^{*}=',num2str(maxrd),'$' newline '$C_r^{*}=' num2str(target_ratio) '$'];  
% draw arrow
annotation('arrow',[ start_point_x,end_point_x],[start_point_y,end_point_y]);           
% add text
annotation('textbox',[ start_point_x-0.13,start_point_y+0.024,0.1,0.1],'String',f2text,'LineStyle','none','Interpreter','latex','FitBoxToText','on'); 
% mark target 
