function [newSAT, savest] = thetaThresholdCircularMaze(obj,spikearraytracked,timeWindow,thetaThresh,plotThetaHist)
%send in methodsfigure for instance, and then the spikearraytracked object.
%see plotPlaceMaps1 in MethodsFigure.
%output is a new spikearraytracked object with only spiketimes (in
%spikearraytrack.SpikeTableInSamples) that are in high theta periods.

if ~exist("thetaThresh","var")
   thetaThresh = 0; %ideally want point that splits a bimodal curve. plot theta hist below to get sense of data.
end

%position already has the time window
lfpWindow = obj.LFP.getTimeWindow(timeWindow);
satWindow = spikearraytracked.getTimeWindow(timeWindow);         

x = obj.positionData.data.X;
z = obj.positionData.data.Z;

posSR = obj.positionData.time.getSampleRate;
lfpSR = lfpWindow.TimeIntervalCombined.getSampleRate;
spikeSR = satWindow.TimeIntervalCombined.getSampleRate;

% posZT = obj.positionData.time.getTimePointsZT;
% lfpZT = lfpWindow.TimeIntervalCombined.getTimePointsZT;
% spikeZT = satWindow.getSpikeTimesZT;
% 
% posAbs = obj.positionData.time.getTimePointsInAbsoluteTimes;
% lfpAbs = lfpWindow.TimeIntervalCombined.getTimePointsInAbsoluteTimes;
% spikeAbs = satWindow.getSpikeTimes;

%DON"T do just samples because that doesn't get precise matching if there
%is more than one file.
%posTimepoints = obj.positionData.time.getTimePointsInSamples/posSR;
%lfpTimepoints = lfpWindow.TimeIntervalCombined.getTimePointsInSamples/lfpSR;
%spiketimes = double(satWindow.SpikeTableInSamples.SpikeTimes)/spikeSR;

posTimepoints = seconds(obj.positionData.time.getTimePointsInAbsoluteTimes-timeWindow(1));
lfpTimepoints = seconds(lfpWindow.TimeIntervalCombined.getTimePointsInAbsoluteTimes - timeWindow(1));
spiketimes = double(seconds(satWindow.getSpikeTimesInAbsoluteTime-timeWindow(1)));

[thetaPow,thetaTfm] = getPowerTimeWindow(obj,timeWindow,[5 12]);
%lfpTime = obj.LFP.TimeIntervalCombined.getTimePointsInSamples/lfpSR;
thetaBinned = interp1(lfpTimepoints,thetaPow.Values,posTimepoints,'pchip'); %find theta power at the times where lfpTime = binnedPos time.
%pchip because lfpTimepoints and posTimepoints don't precisely match
zTheta = zscore(thetaBinned);
% 
% if plotThetaHist
%     thetaHist = histogram(zTheta,100);
%     figure
%    % plot(thetaHist)
% end

%abovethresh = find(zTheta > thetaThresh); 
abovethresh = find(thetaBinned > median(thetaBinned));
abovethreshTimepoints = posTimepoints(abovethresh);

%belowthresh = find(thetaBinned < median(thetaBinned));
%belowthreshTimepoints = posTimepoints(belowthresh);

% Only save spiketimes > threshold

i = 1;
for s = 1:length(spiketimes)
    if any(abs(abovethreshTimepoints-spiketimes(s)) < 1/posSR)
        savest(i) = s; %get index for spikes to save
        i = i + 1;
    end
end

SpikeTimes = satWindow.SpikeTableInSamples.SpikeTimes(savest);
SpikeCluster = satWindow.SpikeTableInSamples.SpikeCluster(savest);
% 
% 
% X = spikearraytracked.TimesInSamples.X(savest);
% Y = spikearraytracked.TimesInSamples.Y(savest);
% Z = spikearraytracked.TimesInSamples.Z(savest);

newSpikeTableInSamples= table(SpikeTimes,SpikeCluster);

newSAT = satWindow;
newSAT.SpikeTableInSamples = newSpikeTableInSamples;