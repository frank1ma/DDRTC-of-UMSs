
%  Flat Output Search & Relative Degree ID - Script for Relatek
%  Date : 07 - 12 - 2019
%  Frank(Shangjie) Ma
%  ----------------------------------------------------------------------
%   - parallel computing config.
%   - Simulation & model selection
%   - Args config.
%   - RelaTek
%   - Plots

% clear 
% clc
%%
% parallel cpu computing
% localpc = parcluster('local');
% parpool('local',localpc.NumWorkers)
%%
%Simulation Config.
% 
% % sampling
% fs = 1000;  % Sample frequency
% tval = 20; % time for simulation
% 
% % reference
% dt = 1 / fs;         % step siz
% t = (0:dt:tval-dt)'; % time vector
% refin =  0.1*randn(1,size(t,1));   % random input
% %load('sim1_refin.mat','refin')
% %refin = chirp(t,1,t(end),500);
% % % linear output matrix
% % C = [1 0 0;0 0 1];
%  C = [1 0 0 0;0 0 1 0];
%%
% Model Selection for simulation

%[out,uin,t] = rdid_lnr_mdl_ivp(fs,tval,C,refin);      % inverted pendulum
%[out,uin,t] = rdid_lnr_mdl_3dsys(fs,tval,C,refin);    % 3d system
%[out,uin,t] = rdid_lnr_mdl_3dsp(fs,tval,C,refin);     % special 3d system
%[out,uin,t] = rdid_nonlnr_mdl_4d(fs,tval,C,refin);     % nonlinear 4d system
%[out,uin,t] = rdid_lnr_mdl_flexible_link(fs,tval,C,refin);
%%
%Exp. data
load('theta_exp_1k_60s_p1_use.mat')
load('alpha_exp_1k_60s_p1_use.mat')
% time vector
 t = theta(1,1:60000)';
%sampling
 fs = 1000;  % Sampling frequency
%reference
 dt = 1 / fs;         % step size

uin = theta(2,1:60000)';  % random input
out = [theta(3,1:60000);alpha(2,1:60000)];
plot(t,out(1,:),t,out(2,:))
%%
% Args Config.

% number of dimensions not required
% it is optional in function getrd to limit the highest relative degree 
 n = 4;

raw_data = out';
raw_input = uin;

% search settings for combination
srch_opts= {[-2,2],1/1000};  % {[start end],search step}
len_srch = size([srch_opts{1}(1):srch_opts{2}:srch_opts{1}(2)],2);

% data interval args
% mode = 'auto' or 'manul'
% for 'auto', a interval should be specified 
% arg list: {{mode,[interval(%)]}}
%crit_opts = {'manul',[]};
crit_opts = {'auto',[0.001 0.06]};

% optional windowing function. Default is Kasiver with beta 35
% win_opts = kaiser(size(t,1)/2+1,35);

win_opts = kaiser(5000,35);

% optional filter applying to raw output and input. Defalut empty.
% see MATLAB Signal Processing Toolbox
filt_opts =[];
%
plot_opts = -0.1;
%%
% core detection function - RelatTek
% returns ratios with relative degree of chosen output combinations
%========================================================================
tic
[rdval,combval,ptgval] = RelaTek(n,t,raw_data,raw_input,srch_opts,crit_opts,win_opts,filt_opts,plot_opts);
%[rdval,combval,ptgval] = RelaTek(n,t,raw_data,raw_input,srch_opts,crit_opts,win_opts,filt_opts);
toc
%%
% Plots
%========================================================================
maxrd = max(rdval);             % max rd
minrd = min(rdval);             
num_of_max_rd = find(rdval==maxrd); % possible ratio for max rd
if size(combval,1) ~= 1         % manul mode
    if (maxrd~=minrd)&&(size(num_of_max_rd,1)>1) 
        
        % confidence rate
        ptgcfdt= ptgval(num_of_max_rd)./sum(ptgval(num_of_max_rd));
        idx_of_max_ptgcfdt = find(ptgcfdt==max(ptgcfdt));
         
        target_list = combval(num_of_max_rd(1) + idx_of_max_ptgcfdt - 1);
        [~,target_idx] = min(abs(target_list));
        target_ratio = target_list(target_idx);
        figure(2)
        yyaxis left
        plot(combval,rdval)
        ylabel('Relative Degree')
        xlabel('Ratio')
        yyaxis right
        bar(combval(num_of_max_rd),ptgval(num_of_max_rd))
        ylabel('Confidence Rate')
        title('Relative Degree and Confidence Rate')
        
        figure(3)
        yyaxis left
        plot(combval(num_of_max_rd(1)-2:num_of_max_rd(end)),rdval(num_of_max_rd(1)-2:num_of_max_rd(end)))
        
        ylabel('Relative Degree')
        xlabel('Ratio')
        yyaxis right
        bar(combval(num_of_max_rd),ptgval(num_of_max_rd))
        ylabel('Confidence Rate')
        title('Relative Degree and Confidence Rate(Zoomed in)')
       
        target_ratio
        maxrd
        
    else                                         % if only one max rd 
        figure
        plot(combval,rdval)
        ylabel('Relative Degree')
        xlabel('Ratio')
        title('Relative Degree')
        
        target_ratio = combval(num_of_max_rd)
        maxrd
        
    end
end