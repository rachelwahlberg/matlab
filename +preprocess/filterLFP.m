function filteredlfp = filterLFP(rawLFP,fs,varargin)
% OLD is a straight box filter. WIll be changing to Gaussian in the new rendition! 

p = inputParser;

addRequired(p,'rawLFP',@isnumeric)
addParameter(p,'sr',1250,@isnumeric)% sampling rate
addParameter(p,'filtervalues',[1,625],@isnumeric)
addParameter(p,'remove60hz',1,@islogical)
addParameters(p,'sumacrosselecs',1,@islogical)

parse(p,varargin{:})
rawLFP = p.Results.rawLFP;
sr = p.Results.sr;
filtervalues = p.Results.filtervalues;
remove60hz = p.Results.remove60hz;
sumacrosselecs = p.Results.sumacrosselecs;

[Npoints, Nelec] = size(rawLFP);
dF = 1000/Npoints; 

rawfreq = zeros(1,ceil(Npoints/2)); % rawfreq =-Fs/2:dF:Fs/2-dF; 
rawfreq = 0:dF:sr/2-dF; %only includes the positive frequencies. 0 to nyquist
  
% f = linspace(0,Fs,Npoints); % useful to check out the filters. 0 to nyq,-nyq to 0.
     
%% Employ functions made below

        %%%%% wide bandpass filter - all frequencies (up to 100 or 500, or others)
        pass.hi = filtervalues(1); %%%% everything below 3, filter DC and other low rhythms
        pass.lo = filtervalues(2);  % if 500, nyq = Npoints/2;
        pass.sixtyhzcut = 1;
        
        filt = getfilter(pass);
        ht = getht(filt);
        
        filteredlfp = real(ht);
%% Make functions for making filter and getting hilbert

    function filt = getfilter(pass)
        filt = zeros(Npoints,1);
        
        fhi=floor(pass.hi/sr *Npoints); % high pass
        flo=ceil(pass.lo/sr *Npoints); % low pass
        higauss = normcdf(1:Npoints,fhi,(2/sr)*Npoints)'; %2/Fs originally in Tay's code. 1/Fs is 1 hz
        logauss = -normcdf(1:Npoints,flo,(2/sr)*Npoints)'; %1/Fs gives sharper cut off but I'm still worried about including some DC

        filt = higauss + logauss; % this does not pass any negative frequencies. will result in hilbert. 
%         filt(fhi:flo) = 1;
        
        if pass.sixtyhzcut == 1;
            hi60 = floor(59.95/sr *Npoints); % filter 60 Hz +/- 0.05 Hz
            lo60 = ceil(60.05/sr *Npoints);

            filt(hi60:lo60) = 0; % don't need to filt neg freqs cause only the positive are being passed 
        end
    end

 function ht = getht(filt)
        
        if sumacrosselecs == 1
            % sum traces first, then run the hilbert
            rawLFPsum = zeros(Npoints,1);% This can ONLY look at the story across the entire array, not ind elecs
            for n = 1:Nelec
                rawLFPsum = rawLFP(:,n)+rawLFPsum;
            end
            
            tmpFFT = fft(rawLFPsum);
            adjFFT = abs(tmpFFT).*exp(1i.*(angle(tmpFFT))); % shift each frequency's phase by -phaseDistort
            
            ht = (ifft(adjFFT.*filt*2))/Nelec; %gives you complex result
            
        elseif sumacrosselecs == 0 % looking at each elec individually, the last input is neuron/elecIDs
            for n = 1:Nelec;
                tmpFFT = fft(rawLFP(:,n));
                adjFFT = abs(tmpFFT).*exp(1i.*(angle(tmpFFT)));
                ht(:,n) = (ifft(adjFFT.*filt*2));
            end
        end
    end
end
