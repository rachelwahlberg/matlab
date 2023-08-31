function [wideht,narrowht] = htfilter(rawLFP,fs,varargin)
% OLD is a straight box filter. WIll be changing to Gaussian in the new rendition! 

[Npoints, Nelec] = size(rawLFP);
Fs = 1000; %sampling freq
dF = 1000/Npoints; 

Narg = nargin;
if Narg == 3
IDs = varargin{1};
end

%%% adjust transfer function - phases at low freqs are distorted due to a Cerebrus issue.
% freq/phase values come directly from Michael Okun via email, based off of Okun, 2017 ("Artefactual origin...")
freq =  [0.0291    0.0601    0.1209    0.1911    0.2546    0.3616    0.8161    1.1161    2.0446    3.0321   11.2444   20.6612   30.3951];
phase = [2.2151    1.8024    1.3725    1.0934    0.9126    0.6852    0.3595    0.2642    0.1891    0.0946    0.0218   -0.0064   -0.0262];
rawfreq = zeros(1,ceil(Npoints/2));
% rawfreq =-Fs/2:dF:Fs/2-dF; 

rawfreq = 0:dF:Fs/2-dF; %only includes the positive frequencies. 0 to nyquist
freqUse = rawfreq >= freq(1) & rawfreq <= freq(end); % frequencies for which phase distortion was measured
phaseShift = zeros(Npoints,1);
phaseShift(freqUse) = interp1(log(freq),phase,log(rawfreq(freqUse)),'pchip');
phaseShift(~freqUse) = 0; % if not in that frequencies, phase shift is 0
  
% f = linspace(0,Fs,Npoints); % useful to check out the filters. 0 to nyq,-nyq to 0.
     
%% Employ functions made below

        %%%%% wide bandpass filter - all frequencies (up to 100 or 500, or others)
        widepass.hi = fs.hipass; %%%% everything below 3, filter DC and other low rhythms
        widepass.lo = fs.lopass;  % if 500, nyq = Npoints/2;
        widepass.sixtyhzcut = 1;
        
        widefilt = getwidefilt(widepass);
        wideht = getht(widefilt);
        
        %%%%% narrow bandpass filter - to get at low frequency oscillations
        narrowpass.hi = fs.narrowhi;
        narrowpass.lo = fs.narrowlo;
        narrowpass.sixtyhzcut = 0;
        
        narrowfilt = getnarrowfilt(narrowpass);
        narrowht = getht(narrowfilt);

%% Make functions for making filter and getting hilbert

    function filt = getnarrowfilt(pass) 
        fhi=floor(pass.hi/Fs *Npoints); % high pass val
        flo=ceil(pass.lo/Fs *Npoints); % low pass val
        fmu = round((fhi + flo)/2);
        
        fsig = (flo - fmu)/2; % the filter bands are 2sd from mu
        filt = normpdf(1:Npoints,fmu,fsig)';
        filt = filt .* (1/max(filt)); % to get the max to be one
    end

    function filt = getwidefilt(pass)
        filt = zeros(Npoints,1);
        
        fhi=floor(pass.hi/Fs *Npoints); % high pass
        flo=ceil(pass.lo/Fs *Npoints); % low pass
        higauss = normcdf(1:Npoints,fhi,(2/Fs)*Npoints)'; %2/Fs originally in Tay's code. 1/Fs is 1 hz
        logauss = -normcdf(1:Npoints,flo,(2/Fs)*Npoints)'; %1/Fs gives sharper cut off but I'm still worried about including some DC

        filt = higauss + logauss; % this does not pass any negative frequencies. will result in hilbert. 
%         filt(fhi:flo) = 1;
        
        if pass.sixtyhzcut == 1;
            hi60 = floor(59.95/Fs *Npoints); % filter 60 Hz +/- 0.05 Hz
            lo60 = ceil(60.05/Fs *Npoints);

            filt(hi60:lo60) = 0; % don't need to filt neg freqs cause only the positive are being passed 
        end
    end
end
