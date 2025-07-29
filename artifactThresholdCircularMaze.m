function [newSAT,savest] = artifactThresholdCircularMaze(obj,spikearraytracked,timeWindow)
%send in methodsfigure for instance, and then the spikearraytracked object.
%see plotPlaceMaps1 in MethodsFigure.
%output is a new spikearraytracked object with only spiketimes (in
%spikearraytrack.SpikeTableInSamples) that are below LFP artifact
%threshold.
%position already has the time window

lfpWindow = obj.LFP.getTimeWindow(timeWindow);
satWindow = spikearraytracked.getTimeWindow(timeWindow);

x = obj.positionData.data.X;
z = obj.positionData.data.Z;

posSR = obj.positionData.time.getSampleRate;
lfpSR = lfpWindow.TimeIntervalCombined.getSampleRate;
spikeSR = satWindow.TimeIntervalCombined.getSampleRate;

posTimepoints = seconds(obj.positionData.time.getTimePointsInAbsoluteTimes-timeWindow(1));
lfpTimepoints = seconds(lfpWindow.TimeIntervalCombined.getTimePointsInAbsoluteTimes - timeWindow(1));
spiketimes = double(seconds(satWindow.getSpikeTimesInAbsoluteTime-timeWindow(1)));

lfpBinned = interp1(lfpTimepoints,lfpWindow.Values,posTimepoints,'pchip');
zLFP = zscore(lfpBinned);

%abovethresh = find(zTheta > thetaThresh); 
belowthresh = find(abs(zLFP) < 5); %above thresh is bad
belowthreshTimepoints = posTimepoints(belowthresh);

abovethresh = find(abs(zLFP) > 5); %above thresh is bad
abovethreshTimepoints = posTimepoints(abovethresh);

% Only save spiketimes > threshold

i = 1;
for s = 1:length(spiketimes)
    if any(abs(belowthreshTimepoints-spiketimes(s)) < 1/posSR)
        savest(i) = s; %get index for spikes to save
        i = i + 1;
    end
end

SpikeTimes = satWindow.SpikeTableInSamples.SpikeTimes(savest);
SpikeCluster = satWindow.SpikeTableInSamples.SpikeCluster(savest);

newSpikeTableInSamples= table(SpikeTimes,SpikeCluster);

newSAT = satWindow;
newSAT.SpikeTableInSamples = newSpikeTableInSamples;