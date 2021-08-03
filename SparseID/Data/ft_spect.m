function [freq,Am,Ph]=ft_spect(dataIn,dt,varargin)
% ft_spect (version 2.0) calculates Amplitude and Phase spectra of an input
% signal with the desired frequency resolution and also filters the Phase 
% spectrum for suppressing the floating rounding-off error.
% 
% NOTICE#1: ft_spect can NOT remove the spectral leakage.
% 
% NOTICE#2: Discrete Fourier transform (DFT) looks at the input signal as one
% period of a periodic signal and discretizes the frequency spectrum of this 
% periodic signal based on the length of the input signal. For a signal with 
% sampling frequency Fs, over the time of T=N?t, the frequency bins (a.k.a 
% frequency resolution in the meaning of distinguishing frequency of f1 from
% f2) are spaced ?f=1/T=Fs/N; thus, the frequency resolution of DFT only depends
% on the length of the input signal (T). BUT, zero-padding does NOT increase 
% the frequency resolution and does NOT reveal more information about the 
% spectrum, it only interpolates amplitudes between bins. For increasing the 
% spectral resolution, a long duration of measurements is necessary, because 
% DFT looks at the input signal as one period of a periodic signal; therefore,
% repeating the input signal is acceptable and doesn't produce any artifacts.
% BUT, in this case, the length of the input signal is increased, and
% consequently, the spectral resolution also increases.
% 
% NOTICE#3: Phase spectrum because of floating rounding-off error is very noisy.
% Small rounding-off error in the "arctan" calculation produces significant noise 
% in the result of the phase spectrum. For suppressing this kind of noise, ft_spect 
% uses a threshold filtering. It means if the amplitude of specific frequency 
% is less than the predefined threshold value, it put zero instead of it.
% 
% IF YOU USE THIS PROGRAM IN YOUR RESEARCH, PLEASE CITE THE FOLLOWING PAPER:
% 
% Afshin Aghayan, Priyank Jaiswal, and Hamid Reza Siahkoohi (2016).
% "Seismic denoising using the redundant lifting scheme." GEOPHYSICS, 81(3), V249-V260.
% https://doi.org/10.1190/geo2015-0601.1
% 
% ** Please share your suggestions and idea for improving ft_spect
% through aghayan@okstate.edu or afshin.aghayan@gmail.com
% 
% Version 1.0 (Spring 2017) ft_spect v1.0 is written and tested in MATLAB R2013a.
% 
% Version 2.0 (Spring 2020) the following changes are applied:
%             1) You can define your desired frequency resolution (?f)
%             2) It is much faster than the v1.0
%             3) Added an example  for comparing the output of the usual FFT
%               and this program (Only type ft_spect for DEMO; look at 
%               ft_spect_demo function at the end of the program for more details)
%              
% Afshin Aghayan
% afshin.aghayan@gmail.com
% 405-334-7184
% 
% --------------------------------INPUT------------------------------------
% 
% By typing ft_spect without any argument, the program runs a DEMO test; 
% look at ft_spect_demo function at the end of the program for more details
% 
% 
% dataIn        The 1D input signal, it can be a row or column vector. Its
%               value must be a real number.
% 
% dt            Time interval or time sampling of the dataIn. It should be
%               a number with the second unit.
% 
% varargin      Combination of OPTIONAL keywords (KeyName) and their values
%               (value_keyName):
%               e.g. []=ft_spect(dataIn,dt,'keyName#1',value_keyName#1, ...)  
%               The first argument is the keyword (it must be a string),
%               and the next argument is its value (could be a string, a
%               number, a vector, or a cell). Each keyword has a default
%               value. 
%               
% POSSIBLE KEYWORDS ARE:
% 
%   'plot'      For defining with spectrum should be display; its value is
%               a string and it has three options:
%               1) 'a' means only display Amplitude spectrum
%               2) 'p' means only display Phase spectrum
%               3) 'ap' means display both amplitude and phase spectra
%               e.g. [freq,Am,Ph]=ft_spect(dataIn,dt,'plot','ap')
%               Default: 'ap'
% 
%   'resolution' For defining frequency resolution (?f). its value must be
%               a number and its unit is Hz.
%               e.g. [freq,Am,Ph]=ft_spect(dataIn,dt,'resolution',0.001)
%               Default: 0.1 Hz
% 
%   'threshold' For defining the threshold value for suppressing floating 
%               rounding-off error noise. All amplitudes that their absolute
%               value are less than this value become zero. Its value must
%               be a number.
%               e.g. [freq,Am,Ph]=ft_spect(dataIn,dt,'threshold',0.1)
%               Default: 0 (means no thresholding)
% 
%   'db'        Display the Amplitude spectrum in decibels. Its value must 
%               be a string and it has two options:
%               1)'y' means convert magnitude to decibels (dB)
%               2)'n' means don't convert magnitude to decibels
%               e.g. [freq,Am,Ph]=ft_spect(dataIn,dt,'db','y')
%               Default: 'n'
% 
%   'deg'       Display the Phase spectrum in degree. Its value must be a 
%               string and it has two options:
%               1)'y' means convert radian to degree
%               2)'n' means show phase based on radian
%               e.g. [freq,Am,Ph]=ft_spect(dataIn,dt,'deg','n')
%               Default: 'y'
% 
%   'figure'    Display the spectrum in a new window. Its value must be a 
%               string and it has two options:
%               1)'y' means create a new window
%               2)'n' means don't create a new window
%               e.g. [freq,Am,Ph]=ft_spect(dataIn,dt,'figure','n')
%               Default: 'y'
% 
%   'color'     Definig the line color in amplitude and phase spectra. Its 
%               value must be a string between the following options:
%               'b' for blue; 'r' for red; 'g' for green; 'k' for black
%               'y' for yellow; 'm' for magenta; 'c' for cyan
%               e.g. [freq,Am,Ph]=ft_spect(dataIn,dt,'color','r')
%               Default: 'b'
% 
%   'title'     Defining a title for amplitude and phase spectra. Its value
%               must be a string.
%               e.g. [freq,Am,Ph]=ft_spect(dataIn,dt,'title','this is a test!')
%               Default: for the amplitude spectrum->'Amplitude spectrum'
%                        for the phase spectrum->'Phase spectrum'
% 
% 
% -----------------------------OUTPUTS-------------------------------------
% 
%   freq        Frequency vector
%   Amp         Amplitude vector
%   Ph          Phase vector

%% --------------------------------DEMO------------------------------------
if nargin<1, ft_spect_demo; return; end

%% --------------------------SET DEFAULT VALUES----------------------------
varName.plot='ap';
varName.resolution=0.1;
varName.threshold=0;
varName.db='n';
varName.deg='y';
varName.figure='y';
varName.color='b';
varName.title=0;

% -------------------------CHECK OPTIONAL INPUTS---------------------------
for i=1:2:length(varargin) 
     if isfield(varName,varargin{i}) %the first argument must be a keyword
         varName.(varargin{i})=varargin{i+1}; % *
     else
        error([' Non-existent keyword: ' varargin{i} ]) 
     end
end

% * the second argument is variable value

%% ------------------------------Main Code---------------------------------

% Make sure dataIn alwaye be a 1D vector and convert row to column vector 
if size(dataIn,1)>=2 && size(dataIn,2)>=2
    error('Input signal should be a 1D signal!')
elseif size(dataIn,1)<2
    dataIn=dataIn';
end

[m,~]=size(dataIn);
Fs =1/dt;               % Sampling frequency

% calculating bin number based on the desire resolution
binNumber=Fs/varName.resolution;
% extending dataIn for reaching to the desired resolution
super_signal=repmat(dataIn,ceil(binNumber/m),1);

%% ---------------------------Fourier Transform----------------------------
 NFFT=length(super_signal);                
 ft =(1/NFFT)*fft(super_signal);                    
 freq = Fs*(0:(NFFT/2))/NFFT;
 Am=abs(ft);
 
%% --------------------------Applying threshold----------------------------
if varName.threshold==0
    Ph=angle(ft);
else
    ft(Am<varName.threshold)=0;
    Ph=angle(ft);
    Am(Am<varName.threshold)=0;
end

if strcmp(varName.deg,'y')
    Ph=radtodeg(Ph);
end

%% -----------------------------Display------------------------------------

% -------------------------Amplitude spectrum------------------------------
if strcmp(varName.plot,'ap') || strcmp(varName.plot, 'a')
    
    if strcmp(varName.figure,'y'); figure; end
    if strcmp(varName.db,'n')
        plot(freq,2*Am(1:NFFT/2+1),varName.color);
        ylabel('Amplitude')
    else
        plot(freq,mag2db(2*Am(1:NFFT/2+1)),color)
        ylabel('Amplitude (dB)')
    end
    xlim([-inf max(freq)])
    xlabel('Frequency (Hz)')
    
    if varName.title~=0
        title(varName.title)
    else
        title('Amplitude spectrum')
    end
end
    
% --------------------------Phase Spectrum---------------------------------
if (strcmp(varName.plot,'ap') || strcmp(varName.plot, 'p'))
    
    if (strcmp(varName.plot,'ap') || strcmp(varName.figure,'y')); figure; end
    plot(freq,Ph(1:NFFT/2+1),varName.color)
    
    if varName.title~=0
        title(varName.title)
    else
        title('Phase spectrum')
    end
    
    xlabel('Frequency (Hz)')
    if strcmp(varName.deg,'y')
        ylabel('Angle (deg)');ylim([-180 180])
    else
        ylabel('Angle (rad)');ylim([-pi pi])
    end
       
end

%% -----------------------------DEMO---------------------------------------
    function ft_spect_demo
        
        clear
        clc
        % --------------------------test signal----------------------------
        dts=0.032;
        L = 125; 
        frequency=[2 4.25 8 12.5 13.25];
        amplitude=[5.6 0.8 0.3 0.75 2.1];
        phase=[-57 30 74 50 -121];
        
        Fsample = 1/dts;            % Sampling frequency
        t = (0:L-1)*dts;
        s=@(i) amplitude(i)*cos(2*pi*frequency(i)*t+(phase(i)*pi/180));
        
        signal=0;
        for j=1:length(frequency)
            signal=signal+s(j);
        end
                   
        % ----------------------------Regular FFT--------------------------
        % NOTICE: Applied as the same as in the "MATLAB 2019b help" 
        f = Fsample*(0:(L/2))/L;
        Y = fft(signal);
        P2 = abs(Y/L);
        P1 = P2(1:floor(L/2)+1);
        P1(2:end-1) = 2*P1(2:end-1);

        % ---------------------Amplitude spectrum--------------------------
        subplot(2,1,1)
        plot(frequency,amplitude,'*')
        hold on,plot(f,P1,'g');xlim([-inf max(f)]) 
        hold on,ft_spect(signal,dts,'plot','a','resolution',0.001,'figure','n','color','r');
        legend('Real value','Regular FFT','ft spect; resolution=0.001 Hz')
        
        % -----------------------Phase spectrum----------------------------
        subplot(2,1,2)
        plot(frequency,phase,'*')
        hold on; plot(f,radtodeg(angle(Y(1:floor(L/2)+1))),'g'); xlim([-inf max(f)])
        hold on,ft_spect(signal,dts,'plot','p','resolution',0.001,'figure','n','color','r','threshold',0.01);
        legend('Real value','Regular FFT','ft spect; threshold value=0.1')

    end

end