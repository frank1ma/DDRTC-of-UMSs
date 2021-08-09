function [max_rd,ptgOut] = getrd(dim,estf,faxis,min_num_seg)
%getrd is a function to compute slope of power spectral density near
% chosen interval. The calculation of slope is affected by interval and
% outlier method as well as the threshold factor.
%
%   dim       - preserved for future use
%   psd_data  - abs of power spectral density of outputs
%   psd_u     - abs of power spectral density of inputs
%   f         - frequncies of psd
%   coeff     - coefficients of combination
%   crit_opts - criteria for slope
%
%   Frank S. Ma

% Max Relative Degree
% Covering Index
max_rd = 0;
ptgOut = 0;

% Uncomment to check estf within interval
% plot(faxis,estf)

% divide the estf into portions in terms of peaks and local minima
[seg]=divpeak(estf,faxis);

for k = 1:size(seg,2)
    if seg(k).ptg < min_num_seg
        continue;
    end
    val = polyfit(seg(k).faxis,seg(k).estf,1);
    slope_of_seg = abs(val(1));
    rd_candiate =round(slope_of_seg);
    if rd_candiate > max_rd
        max_rd = rd_candiate;
        ptgOut = seg(k).ptg;%seg(k).ptg;
    end
end

% if max_rd > dim    % max rd shoud not exceed dim.
%     max_rd = dim;
% end

end

%[EOF]