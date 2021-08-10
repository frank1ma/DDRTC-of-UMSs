function [max_rd,ptgOut] = getrd(dim,estf,faxis,min_num_seg,plot_opts,coeff)
%getrd is a function to compute slope of power spectral density within
% chosen interval. The calculation of slope is affected by selection of interval 
% and outlier method, as well as the threshold factor.
%
%   dim      - preserved for future use, dimension information                                   
%   estf      - data points of estimated frequency response of trial output
%                 within chosen frequency band
%   faxis     - frequency axis
%   min_num_seg  - threshold for valid subset of estf.
%
%   Frank S. Ma
%   Last modification :  10-21-2019

% Max Relative Degree
% Main Result Variable Init
max_rd = 0;
ptgOut = 0;
max_seg = 0;
% Uncomment to check estf within interval
% plot(faxis,estf)
if(plot_opts==coeff)
    plot_flag = 1;
else
    plot_flag = 0;
end
 
% divide the estf into serveral subsets in terms of peaks and local minima
[seg]=divpeak(estf,faxis,plot_flag);

% filter ineligible subset out of the candidates
for k = 1:size(seg,2)
    
    if seg(k).ptg < min_num_seg   % number of data points under threshold 
        continue;
    end
    
    val = polyfit(seg(k).faxis,seg(k).estf,1);  % estimate slope
    
    if val(1) > 0                                            % if the slope is positive
        continue;
    end
    
    seg(k).slope_of_seg = abs(val(1));                    % take absolute value
    
    rd_candidate =round(seg(k).slope_of_seg);    % round off slope to integer
    
    if rd_candidate > max_rd                      % take largest relative degree and its size of subset.
        max_rd = rd_candidate;
        %ptgOut = abs((rd_candidate - slope_of_seg)/rd_candidate); 
        ptgOut = seg(k).ptg;
        max_seg = k;
    end
end

% if max_rd > dim    % max rd shoud not exceed dim.
%     max_rd = dim;
% end
if(plot_flag)
save('segout.mat')
end
% uncomment to check maximum relative degree
% max_rd;                  
end

%[EOF]