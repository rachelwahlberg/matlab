function [rawdata,tmpbad] = glitchdetector(varargin)

% Removes dropout and artifact periods
% Generates file of [start,end] bad timeperiods (removetimestamps)
%% Inputs
% Default values
p = inputParser;
addRequired(p,'data',@isnumeric) %or ismatrix?
addRequired(p,'lfpsamplingrate',@isnumeric)
addParameter(p,'winsize',100,@isnumeric)
addParameter(p,'minInterval',5000,@isnumeric)
addParameter(p,'saveMat',false,@islogical);
addParameter(p,'percentiles',[2.2 99.5],@isnumeric) % defaults manually chosen
addParameter(p,'eventfilename',[],@ischar)
addParameter(p,'newlfpfilename',[],@ischar)

% assign parameters (either defaults or given)
parse(p,varargin{:})
data = p.Results.data;
lfpsamplingrate = p.Results.lfpsamplingrate;
winsize = p.Results.winsize; %in ms
minInterval = p.Results.minInterval;
saveMat = p.Results.saveMat;
percentiles = p.Results.percentiles;
eventfilename = p.Results.eventfilename;
newlfpfilename = p.Results.newlfpfilename;

%% Create windows, take variances

rawdata = data;
if size(data,2) > 1 % if sent in with multiple channels  
    data = mean(data,2);
end

nWins = floor(length(data)/winsize);
roundTrace= data(1:nWins*winsize);                                            % multiply total number of windows by window size to get the
% full length of the trace, rounded so no leftovers.
shapedTrace=reshape(roundTrace,winsize,length(roundTrace)/winsize);         %gives matrix of win x total number of windows
winVar=var(shapedTrace); % takes variance of each window
pvar = prctile(winVar,percentiles);                                         % need variance rather than straight values because of flat dropouts.
% not using z score because the crazy outliers are pulling the mean up

%% Remove dropout (FLAT parts of signal)

dropoutInd = find(winVar(:)<pvar(1)); %detecting which windows are part of a dropout, if any
start1 = dropoutInd*winsize-99; % first passes - will merged windows.
end1 = dropoutInd*winsize;

if dropoutInd(end) == nWins% fix for the case of it's the last time bin
    end1(end,2) = length(data);
end

i = 1;
drop.start(1) = start1(1);  drop.end = [];
for win = 2:length(dropoutInd)

    % special case of last index
    if win == length(dropoutInd)
        drop.end(i) = end1(win);
        break
    end

    % merge drop periods that are close together
    if start1(win)-end1(win-1)<minInterval % if the dif between win end to next win start is < 500 ms,
        continue
    else
        drop.end(i) = end1(win-1);
        drop.start(i+1) = start1(win);
        i = i + 1;
    end
end

%% Detect artifacts (far ABOVE average activity)

artInd = find(winVar(:)>pvar(2)); %detecting which windows are part of a dropout, if any
start1 = artInd*winsize-99; % first passes - will merged windows
end1 = artInd*winsize;

if artInd(end) == nWins% fix for the case of it's the last time bin
    end1(end) = length(data);
end

i = 1;
art.start(1) = start1(1);  art.end = [];
for win = 2:length(artInd)

    % special case of last index
    if win == length(artInd)
        art.end(i) = end1(win);
        break
    end

    % merge drop periods that are close together
    if start1(win)-end1(win-1)<5000 % if the dif between win end to next win start is < 500 ms,
        continue
    else
        art.end(i) = end1(win-1);
        art.start(i+1) = start1(win);
        i = i + 1;
    end
end

%% Nan output option
% bringing back all channels here
del.start = [drop.start art.start];
del.end = [drop.end art.end];

for d = 1:length(del.start)
    rawdata(del.start(d):del.end(d),:) = nan;
end

%% Timestamp output option in neuroscope format
% has to end in .[3 char identifier].evt
% suggested name, [basename '_removeTimestamps.bad.evt']

tmpbad = [drop.start art.start;drop.end art.end]'; % have start times of bad time periods in col 1, end times in col 2
tmpbad = sort(tmpbad,1,"ascend");

tmpbad2 = tmpbad/lfpsamplingrate * 1000;

events.time = single(vertcat(tmpbad2(:)));

events.description = cell(length(events.time),1);
for i = 1:2:length(events.description)
    events.description{i} = "Start";
    events.description{i+1} = "Stop";
end


if eventfilename % have filename end in .evt
%    WriteEventFileFromTwoColumnEvents (tmpbad2,eventfilename) % Brendon Watson version

   SaveEvents_RW(eventfilename,events) %a buzcode function
end

%% save clean LFP file

if newlfpfilename
    cleanlfp.data = rawdata; % cleaned in "Nan output option" section
    cleanlfp.params.cleanfunction = 'glitchdetector.m';
    cleanlfp.params.percentiles = pvar;
    cleanlfp.params.removetimestamps = del;

    save(newlfpfilename,'cleanlfp')
end


