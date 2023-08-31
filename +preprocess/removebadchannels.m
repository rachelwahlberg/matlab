function [lfpout,badCh] = removebadchannels(lfp)

% takes in STRUCTURE which has .data (nsamps x nchan) and .channels (vec of
% all channels)

% ncan't take in precleaned data if too noisy; 
% can't take in "cleanlfp" as is because of nans
% removing timestamps with any nans 

data = lfp.data;
data(isnan(data(:,1)),:) = []; 

avgCh = mean(double(data),2);

for i=1:size(data,2)
    c = corrcoef(avgCh,double(data(:,i)));
    cc(i) = c(2);
    bad(i) = cc(i)<.7;
end

lfpout = lfp;
lfpout.data(:,bad) = [];
lfpout.channels(bad) = [];

badCh = lfp.channels(bad);
% ch 122 looks like it should be bad too?


