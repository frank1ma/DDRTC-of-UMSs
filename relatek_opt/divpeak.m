function [segout,seg_del,seg_zp,is_line_val]=divpeak(estf,faxis,plot_flag)
% DIVPEAK is a function that divides ESTF into several lines,circumventing
% the local minima and maxima of input data ESTF. DIVPEAK returns several
% segments separated by minima and maxima. Elements of number +- 'range' around
% minma and maxima will be discarded. The subset and size of each subset
% will be returned along with the array of structs.
%
% If no peaks or minima are found, DIVPEAK returns modified estf with left
% side cut off ,where cut-off elements has a quantity of number specified.
%
% Input:
% estf  - estimate transfer function
% faxis - frequency information of estf
%
% Output:
% seg  - array of struct containing the information of segments of estf.
%
% Available for MATLAB CODER
% DATE FOR CREATION : 07-03-2019
% LAST MODIFICATION : 10-21-2019
% FRANK S.MA


% length of estf
% sequence of tf
% left-cut threshold
% logical array for selection
len_estf = size(estf,1);
seq = 1:len_estf;
num_of_left_cut = round(0.05 * len_estf);
TF = true(len_estf,1);
is_line_threshold = 5;

% trend line
p_trend = polyfit([faxis(1);faxis(end)],[estf(1);estf(end)],1);
estf_detrend = estf - polyval(p_trend,faxis);
[maxval,maxidx] = max(abs(estf_detrend));
%is_line_val = abs(maxval/(faxis(end)-faxis(1)));
[p,S] = polyfit(faxis,estf_detrend,1);
is_line_val = S.normr;
%plot(faxis,estf,'g',faxis,estf_detrend,'b',faxis(1),estf(1),'rx',faxis(end),estf(end),'rx',faxis(maxidx),maxval,'o');
%line([faxis(1) faxis(end)],[faxis(1)*p_trend(1)+p_trend(2) faxis(end)*p_trend(1)+p_trend(2)]);
%findpeaks(-estf_detrend,'MinPeakProminence',0.25,'NPeaks',10);

if abs(is_line_val(1)) < is_line_threshold
    TF(1:num_of_left_cut) = false;
    locs_zp=[];
else
    TF(maxidx) = false;
    % find local minima produced by zeros of estf.
    %[~,locs_zeros]=findpeaks(-estf,'MinPeakProminence',0.25,'NPeaks',10);
    %[~,locs_zeros]=findpeaks(-estf_detrend,'MinPeakProminence',0.25,'NPeaks',10);
    [~,locs_zp]=findpeaks((estf),'MinPeakProminence',0.5,'NPeaks',20);
    [~,locs_zp1]=findpeaks(abs(estf),'MinPeakProminence',0.5,'NPeaks',20);
    % find local peak produced by poles of estf.
    %[~,locs_poles]=findpeaks(estf,'MinPeakProminence',0.25,'NPeaks',10);
    %[~,locs_poles]=findpeaks(estf_detrend,'MinPeakProminence',0.25,'NPeaks',10);
    locs_zp = [locs_zp;maxidx;locs_zp1]; 
    % if any local minimum or peak
    %if (~isempty(locs_zeros) || ~isempty(locs_poles)) && (size(locs_zeros,1) < 10 && size(locs_poles,1) < 10)
    if (~isempty(locs_zp)) && (size(locs_zp,1) < 20)
        %locs_zp = [locs_zeros;locs_poles];
        TF(locs_zp)= false;      % mark minima or peak
        %plot(faxis,estf,faxis(~TF),estf(~TF),'kx')
        for i=1:size(locs_zp,1) % mark the neighborhood of minima or peak as false
            loc_val = faxis(locs_zp(i));
            left_bound = max(find(faxis<(loc_val-0.05)));
            right_bound = min(find(faxis>(loc_val+0.05)));
%             cut_off_range_right = round(locs_zp(i)*0.50);            % cutoff range
%             cut_off_range_left = round(locs_zp(i)*0.50);
            %if (locs_zp(i)+cut_off_range_right) <= len_estf
            if (right_bound) <= len_estf
                %TF((locs_zp(i)-cut_off_range_left):(locs_zp(i)+cut_off_range_right))= false; % right not reaches bound
                TF(left_bound:right_bound) = false;
            else
                %TF((locs_zp(i)-cut_off_range_left):len_estf)= false; % right reaches bound
                TF(left_bound:len_estf)= false;
            end
        end
    else
        % no minima and maxima found
        TF(1:num_of_left_cut) = false;
    end
end

% Uncomment to plot region. Red marks are points preserved, blues are
% discarded
% plot(faxis,estf,faxis(TF),estf(TF),'rx')
% plot(faxis,estf,faxis(~TF),estf(~TF),'bx')

% test codes for separation
% TF = logical([1 1 0 0 0 0 0 0 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 1 1 1 1 1]);

% % Uncomment to plot the original estimate transfer function
o_estf = estf;
o_faxis= faxis;
if(plot_flag)
seg_zp.x=o_faxis(locs_zp);
seg_zp.y=o_estf(locs_zp);
seg_del.x=faxis(~TF);
seg_del.y=estf(~TF);
save('seg_del_zp.mat')
end
seg = struct;
estf = estf(TF);             % take points selected
faxis = faxis(TF);

% default output is for only one segment,meaning all points will be taken
% to evaluate the slope. In another words, no obvious peak or local minima.
seg.estf = estf;
seg.faxis = faxis;
seg.ptg = (len_estf - num_of_left_cut)/len_estf;
segout =seg;
% Uncomment to plot tf with  one segment, aka no peaks and local minimum
% plot(o_faxis,o_estf,seg(1).faxis,seg(1).estf,'x')

% diff(seq(NTF))-1, is the position where the number become nonconsecutive
% this is based on the fact that nonconsecutive numbers will show
% difference more than 1.
[~,loc_y]=find((diff(seq(TF))-1)>0);  % showing the border of each segments

% Modify default output if more than one segment needed to evaluate
if ~isempty(loc_y)
    num_of_seg = 1 + size(loc_y,2);
    segarr = repmat(seg,1,num_of_seg);    % creat array of structs
    loc_yy = [0,loc_y,size(estf,1)];   % create indices for separation
    
    for j = 1:num_of_seg
        % get info of each segment
        segarr(j).ptg = (loc_yy(j+1)-loc_yy(j))/ len_estf;
        segarr(j).estf= estf(loc_yy(j)+1:loc_yy(j+1));
        segarr(j).faxis=faxis(loc_yy(j)+1:loc_yy(j+1));
    end
    % Uncomment to plot the intereted segments, example for two segments
    %plot(o_faxis,o_estf,'g-',segarr(1).faxis,segarr(1).estf,'bo',segarr(2).faxis,segarr(2).estf,'rx')
    segout = segarr;
end
end

%[EOF]

