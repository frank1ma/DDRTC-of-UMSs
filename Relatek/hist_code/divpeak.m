function [segout]=divpeak(estf,faxis)
% DIVPEAK is a function that divides EST into several lines,circumventing
% the local minima and maxima of input data EST. DIVPEAK returns several
% segments separated by minima and maxima. Elements of number +- 'range' around
% minma and maxima will be discarded. The percentage of original EST that
% each segment takes up will be returned as ptg.
%
% If no peaks or minima are found, DIVPEAK returns modified estf with left
% side cut off ,where cut-off elements has number specified.
%
% Input:
% estf  - estimate transfer function
% faxis - frequency in terms of estf.
%
% Output:
% seg  - array of struct containing the information of segments of estf.
%
% Available for MATLAB CODER
% DATE : 07-03-2019
% FRANK S.MA


% length of estf
% sequence of tf
% left-cut threshold
% logical array for selection
len_estf = size(estf,1);
seq = 1:len_estf;
num_of_left_cut = round(0.05 * len_estf);
TF = true(len_estf,1);

% find local minima produced by zeros of estf.
[~,locs_zeros]=findpeaks(-estf,'MinPeakProminence',1);
% find local peak produced by poles of estf.
[~,locs_poles]=findpeaks(estf,'MinPeakProminence',1);


% if any local minimum or peak
if ~isempty(locs_zeros) || ~isempty(locs_poles)
    locs_zp = [locs_zeros;locs_poles];
    TF(locs_zp)= false;      % mark minima or peak
    % plot(faxis,estf,faxis(~TF),estf(~TF),'rx')
    for i=1:size(locs_zp,1) % mark the neighborhood of minima or peak as false
        cut_off_range = round(locs_zp(i)*0.35);            % cutoff range
        TF((locs_zp(i)-cut_off_range):(locs_zp(i)+cut_off_range))= false; % left not reaches bound
    end
else
    % no minima and maxima found
    TF(1:num_of_left_cut) = false;
end

% Uncomment to plot region. Red marks are points preserved, blues are
% discarded
% plot(faxis,estf,faxis(TF),estf(TF),'rx')
% plot(faxis,estf,faxis(~TF),estf(~TF),'bx')

% test codes for separation
% TF = logical([1 1 0 0 0 0 0 0 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 1 1 1 1 1]);

% Uncomment to plot the original estimate transfer function
% o_estf = estf;
% o_faxis= faxis;



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
    % plot(o_faxis,o_estf,'g-',seg(1).faxis,seg(1).estf,'bo',seg(2).faxis,seg(2).estf,'rx')
    segout = segarr;
end

end

%[EOF]

