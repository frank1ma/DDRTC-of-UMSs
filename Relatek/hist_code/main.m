%*      *       *       *       *       *       *       *       *       *
%  Script for RelaTeK
%  Date : 07 - 01 - 2019
%  Frank(Shangjie) Ma
%  ----------------------------------------------------------------------
%   - Model Selection
%   - Arguments
%   - RelaTek
%   - Plots
%*      *       *       *       *       *       *       *       *       *

% parallel computing
% localpc = parcluster('local');
% parpool('local',localpc.NumWorkers)

% Model Selection
%========================================================================
RDID_test
%models to be developed :
%test_linear_model
%test_linear_model_inverted_pendulum
%test_nonlinear_model
%test_duffing_model

% Arguments 
%========================================================================
% number of dimensions not required
% it is optional in function getrd to limit the highest relative degree 
% n = [];

t = t';
raw_data = [x1' x3'];
raw_input = u_ref';

% search settings for combination
srch_opts= {[-0.90 -0.89],1/10000};  % {[start end],search step}
len_srch = size([srch_opts{1}(1):srch_opts{2}:srch_opts{1}(2)],2);

% data interval
% mode = 'auto' or 'manul'
% for 'auto', a interval should be specified 
% arg list: {{mode,[interval(%)]}}
%crit_opts = {'manul',[]};

crit_opts = {'auto',[0.01 0.5]};

% optional windowing function. Default is Kasiver with beta 35
% win_opts = kaiser(size(t,1)/2+1,35);

win_opts = [];

% optional filter applying to raw output and input. Defalut empty.
% see MATLAB Signal Processing Toolbox
filt_opts =[];

% core detection function - RelatTek
% returns ratios with relative degree of chosen output combinations
%========================================================================
tic
[rdval,combval,ptgval] = RelaTek([],t,raw_data,raw_input,srch_opts,crit_opts,[],[]);
toc

% Plots
%========================================================================
maxrd = max(rdval);             % max rd
num_of_max_rd = find(rdval==maxrd); % possible ratio for max rd

if size(combval,1) ~= 1
    if size(num_of_max_rd,1) > 1
        
        % confidence rate
        ptgcfdt= ptgval(num_of_max_rd)./sum(ptgval(num_of_max_rd));
        
        figure
        yyaxis left
        plot(combval,rdval)
        ylabel('Relative Degree')
        xlabel('Ratio')
        
        yyaxis right
        bar(combval(num_of_max_rd),ptgcfdt)
        ylabel('Confidence Rate')
        title('Relative Degree and Confidence Rate')
    else
        
        figure
        plot(combval,rdval)
        ylabel('Relative Degree')
        xlabel('Ratio')
        title('Relative Degree')
    end
end