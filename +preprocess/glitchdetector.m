function [cleanlfp,removetimestamps] = glitchdetector(varargin)

% Removes dropout and artifact periods
% Generates file of [start,end] bad timeperiods (removetimestamps)

%%% Outputs:
% cleanlfp is in format of 
%% Inputs
% Default values
p = inputParser;
addParameter(p,'winsec',.3,@isnumeric) % in SECONDS
addParameter(p,'minInterval',5000,@isnumeric)
addParameter(p,'saveMat',false,@islogical);
addParameter(p,'percentiles',[2.2 99.5],@isnumeric) % defaults manually chosen for 11/23
addParameter(p,'eventfilename',[],@ischar)
addParameter(p,'newlfpfilename',[],@ischar)
addParameter(p,'mergelatency',.3,@isnumeric) %300 ms - if wins are less than this val, they are merged

if isnumeric(varargin{1}) % matrix sent in as M samples x N trials
addRequired(p,'data',@isnumeric) %or ismatrix?
addRequired(p,'timestamps',@isnumeric)
addRequired(p,'lfpsamplingrate',@isnumeric)
parse(p,varargin{:})
data = p.Results.data;
timestamps = p.Results.timestamps;
elseif isstruct(varargin{1})
addRequired(p,'lfp',@isstruct) % in as the full lfp structure
addRequired(p,'lfpsamplingrate',@isnumeric)
parse(p,varargin{:})
lfp = p.Results.lfp;
data = lfp.data;
timestamps = lfp.timestamps;

else
    disp('first argument must be either lfp structure of matrix of samples x trials')
end

lfpsamplingrate = p.Results.lfpsamplingrate;
winsec = p.Results.winsec; %in seconds 
minInterval = p.Results.minInterval;
saveMat = p.Results.saveMat;
percentiles = p.Results.percentiles;
eventfilename = p.Results.eventfilename;
newlfpfilename = p.Results.newlfpfilename;
mergelatency = p.Results.mergelatency;
% assign parameters (either defaults or given)

%% Create windows, take variances

rawdata = data;
if size(data,2) > 1 % if sent in with multiple channels 
    clear data
    data = single(mean(rawdata,2));
end

winsize = winsec * lfpsamplingrate; %window in seconds * lfpsamplingrate

nWins = floor(length(data)/winsize);
roundTrace= (data(1:nWins*winsize)).^2;  % multiply total number of windows by window size to get the
% full length of the trace, rounded so no leftovers.
shapedTrace=reshape(roundTrace,winsize,length(roundTrace)/winsize);         %gives matrix of win x total number of windows
winPow = abs(hilbert(single(mean(shapedTrace))));
winVar=var(single(shapedTrace)); % takes variance of each window
pvar = prctile(winVar,percentiles);                                         % need variance rather than straight values because of flat dropouts.
% not using z score because the crazy outliers are pulling the mean up

winPowZ = abs(zscore(winPow));
%winZ = abs(zscore(mean(shapedTrace)));
pZ = 7;
%% Remove flat parts of signal and extra high parts of signal

% first pass - get index of window numbers, convert into start and end
% timepoints (in samples)
%dropoutInd = find(winVar(:)<pvar(1)); %detecting which windows are part of a dropout, if any
%artInd = find(winVar(:)>pvar(2)); 
%firstpass = sort(vertcat(dropoutInd,artInd));

dropoutInd = find(winPowZ>pZ);
%dropoutInd = find(winVar(:)>pvar);
firstpass = dropoutInd;

start1 = firstpass*winsize-(winsize-1); %convert into SAMPLES at starts and ends of windows
end1 = firstpass*winsize;

if firstpass(end) == nWins% fix for the case of it's the last time bin
    end1(end) = length(data);
end

%second pass - merge windows 

i = 1;
secondpass(i,1) = start1(i);
for win = 2:length(firstpass)

    % special case of last index
    if win == length(firstpass)
        secondpass(i,2) = end1(win);
        break
    end

    % merge periods that are close together
    if start1(win)-end1(win-1)<lfpsamplingrate/(1/mergelatency) % if the dif between win end to next win start is < .3 s,
        continue
    else
        secondpass(i,2) = end1(win-1);
        secondpass(i+1,1) = start1(win);
        i = i + 1;
    end
end

% third pass - deleting too small of windows, extending the larger windows
i = 1;
for ep = 1:length(secondpass)
    if secondpass(ep,2) - secondpass(ep,1) < lfpsamplingrate/(1/.05) % removing periods that are too short (< 50 ms)
        continue
    else
        thirdpass(i,1) = secondpass(ep,1)-ceil(lfpsamplingrate/(1/0.05)); % extending the longer periods to capture full thing. currently specific to 20230614's artifacts!
        thirdpass(i,2) = secondpass(ep,2)+ceil(lfpsamplingrate/(1/0.2));
        i = i + 1;
    end
end

%% Nan output option
% bringing back all channels here
rawdata = single(rawdata); % otherwise nans become zeros when int16
for d = 1:length(thirdpass)
    rawdata(thirdpass(d,1):thirdpass(d,2),:) = nan;
    timestamps(thirdpass(d,1):thirdpass(d,2),:) = nan;
end

%% Timestamp output option in neuroscope format
% has to end in .[3 char identifier].evt
% suggested name, [basename '_removeTimestamps.bad.evt']

finalpass = thirdpass'/lfpsamplingrate * 1000; % into ms
%finalpassDead = thirdpass'/lfpsamplingrate * 30000;

events.time = single(vertcat(finalpass(:)));

events.description = cell(length(events.time),1);
for i = 1:2:length(events.description)
    events.description{i} = "Start";
    events.description{i+1} = "Stop";
end

if eventfilename % have filename end in .evt
%    WriteEventFileFromTwoColumnEvents (tmpbad2,eventfilename) % Brendon Watson version
   utilities.SaveEvents_RW(eventfilename,events) %a buzcode function
end

if deadfilename %if for PREPHY
    utilities.SaveDeadFile(deadfilename,finalpass) %need to write
end

%% save clean LFP file

removetimestamps.evtformat = finalpass'; % ms from start of LFP file
removetimestamps.sampformat = thirdpass; % in samples, to index into lfptimestamps.fromlfpclockstart
removetimestamps.sampfreq = lfpsamplingrate;

if isnumeric(varargin{1})
    cleanlfp.data = rawdata; % cleaned in "Nan output option" section
    cleanlfp.params.cleanfunction = 'glitchdetector.m';
    cleanlfp.params.percentiles = pvar;
    cleanlfp.params.removetimestamps = finalpass;
    cleanlfp.timestamps = timestamps;
elseif isstruct(varargin{1}) % takes original lfp structure, replaces data with clean data, adds new params
    cleanlfp = lfp;
    cleanlfp.data = rawdata; % cleaned in "Nan output option" section
    cleanlfp.params.cleanfunction = 'glitchdetector.m';
    cleanlfp.params.percentiles = pvar;
    cleanlfp.params.removetimestamps = finalpass;
    cleanlfp.cleantimestamps = timestamps;
    cleanlfp.Filename = [cleanlfp.Filename(1:end-3) 'clean' cleanlfp.Filename(end-3:end)];
    if newlfpfilename
        cleanlfp.Filename = newlfpfilename; %overwrite the filename
        save(newlfpfilename,'cleanlfp')
    end
end




