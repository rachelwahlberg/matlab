%% Preliminary data pipeline
%% 0a) merge files
% Run only pre-sorting! 
usedate = 20240625;
basefolders = preprocess.getbasefolders(usedate); % adapt to include the variables below to not have to rerun each time
mergedfilename = '/data/ExperimentsRSW/CircularMaze/James/20240625/merged_20240625';
%mergedfilename = '/data/ExperimentsRSW/20230614/merged_20230614';
preprocess.OpenEphysPreprocess_RW(basefolders,mergedfilename);

%%% 0b) spyking circus in merged folder - need .params, .prb, .dat, dead.txt (if bad periods)
%%% 0c) phy in merged folder

%% 1a) Define the basepath of the dataset to run. 
% The dataset should at minimum consist of the raw data and spike sorted data.
basepath = '/data/ExperimentsRSW/CircularMaze/20230614/merged_20230614';
basename = 'merged_20230614';
behaviorpath = fullfile(basepath, 'behaviorfiles/');
phyoutputpath =  fullfile(basepath,'home/wahlberg/merged_M1_20211124crs.GUI');
TimeIntervalCombinedpath = fullfile(basepath, [basename '_1250Hz.TimeIntervalCombined.csv']);

cd(basepath)

%% Get neuroscope timepoints
clustering_path = '/home/wahlberg/merged_20230614crs.GUI';
Phy2Neurosuite(basepath,clustering_path)

%% 1b) CE Generate session metadata struct using the template function and display the meta data in a gui

if isfile(fullfile(basepath, [basename '.session.mat']))
    disp('Loading existing session file')
    load([basename '.session.mat'])
else
    session = sessionTemplate(pwd,'showGUI',true); % CELL EXPLORER
    validateSessionStruct(session); % validate the required and optional fields
end

if isfile(fullfile(basepath, [basename '.sessionInfo.mat']))
    disp('Loading existing sessionInfo file')
    load([basename '.sessionInfo.mat'])
else
    sessionInfo = bz_getSessionInfo(basepath); % BUZCODE see < https://github.com/buzsakilab/buzcode/wiki/Data-Formatting-Standards > .
end

%% 2a)  BZ Get the lfp data from the dat file

if isfile(fullfile(basepath, [basename '.lfp']))
    disp('Loading existing LFP file')
    lfp = bz_GetLFP('all','basepath',basepath,'basename',basename);
else
    disp('LFP file does not yet exist, creating one from memmap')
    lfp = bz_GetLFP('all','basepath',basepath,'basename',basename); 
end

%% 2b) BZ/RW Import behavioral tracking from optitrack

if isfile(fullfile(basepath, [basename '.dat.optitrack.behavior.mat']))
    disp('Loading existing optitrack file')
    load(fullfile(basepath, [basename '.dat.optitrack.behavior.mat']));
else
    optitrack = getbehavior(basepath,basename,behaviorpath);
    optitrack = outofmaze(optitrack,'basepath',basepath,'basename',basename);
end

%% 2c) get clean lfp

[cleanlfp,removetimestamps] = glitchdetector(lfp,session.extracellular.srLfp,...
    'percentiles',[0.9 99.93]); %,'saveMat',true,'eventfilename',fullfile(basepath,[basename '.bad.evt']));

%% 2d) Get timestamps of LFP and behavior files relative to start of LFP
% periods with "removetimestamps" will be nans if included

if isfile(fullfile(basepath, [basename '.lfptimestamps.mat'])) && isfile(fullfile(basepath, [basename '.behaviortimestamps.mat']))
    load(fullfile(basepath, [basename '.lfptimestamps.mat']));
    load(fullfile(basepath, [basename '.behaviortimestamps.mat']));
else
    [lfptimestamps,behaviortimestamps] = getTimestamps(basepath,behaviorpath,...
        'TimeIntervalCombinedpath',TimeIntervalCombinedpath,'removetimestamps',removetimestamps);
end

%%%% 
if length(unique(behaviortimestamps.fromlfpfilestart)) < ...
    length(behaviortimestamps.fromlfpfilestart)
    idx = findNonUniqueValues(behaviortimestamps);
    save([basename '_repeatedBehaviorVals_created' date '.mat'],'idx')
    warning('Check your recording for repeated values')
end
    
%startofbehavior = behtimestamps.reltime_merged.Time(1); % startofbehavior relative to start of LFP, in ms. ADD this value
behaviortimestamps.fromlfpfilestart(34200:34400) = ...
    linspace(behaviortimestamps.fromlfpfilestart(34200),...
    behaviortimestamps.fromlfpfilestart(34400),201);

%%%%%%%%%%%%%%%%%% 3) Divide into trials %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 3a) Divide into trials
if isfile([basename '.trials.mat'])
    load([basename '.trials.mat'])
else
    newtrialfilename = 'basename.trials.mat';
trials = gettrials(optitrack,behaviortimestamps,'alternationYmaze',...
    'savetrials',false,'plottrials',true,'basepath',basepath,'basename',basename);
% in clock ms from start of LFP file
end

%% 3b) separate into pre,maze,post epochs 
% focus on the maze part for now

mazeStart = trials.startTrial_clock(1); % in what? 
mazeEnd = trials.endTrial_clock(end);

behstartIndex = find(behaviortimestamps.fromlfpclockstart(:) == mazeStart);
behendIndex = find(behaviortimestamps.fromlfpclockstart(:) == mazeEnd);

lfpstartIndex = behaviortimestamps.indexintolfp(behstartIndex);
lfpendIndex = behaviortimestamps.indexintolfp(behendIndex);

mazelfp = lfp.data(lfpstartIndex:lfpendIndex,:); % Not the cleaned version
premazelfp = lfp.data(1:lfpstartIndex-1,:);
postmazelfp = lfp.data(lfpendIndex+1:end,:);

%% 4 Determine theta power during the active behavior periods
% from the beginning of the first behavior timestamps to the last behavior
% timestamp

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Check EEG states %%%%%%%%%%%%%%%%%%%%%%%%%

testlfp = mazelfp(1000000:end,:); % just pulling out the clean parts for testing
testpre = premazelfp(1:4800000,126); % for same length as testlfp
testpost = postmazelfp(4000000:8800000,126);

CheckEegStates_RW(basepath,basename,'lfp',testpost,'redo_flag',true); % and testpre, and testpost

%%%%%%%%%%%%%Check EEG states with sliding window %%%%%%%%%%%%%%%%%%%%%%%%%
% Put "100" as a window size into WhitenSignal, didn't change much
% figure out what the measurement is for the window

%%%%%%%%%%%%%Check EEG states for entire behavioral period %%%%%%%%%%%%%%%%
CheckEegStates_RW(basepath,basename,'lfp',testlfp,'redo_flag',true); % and testpre, and testpost

%%%%%%%%%%%%%% Try SleepScoreMaster for theta detection %%%%%%%%%%%%%%%%%%

%badCh = session.channelTags.Bad.channels;
%thetaCh = session.channelTags.ThetaChan.channels; % check this

badCh = [0 122];
thetaCh = 126;
SleepState = SleepScoreMaster(basepath,'ignoretime',removetimestamps.sampformat,'rejectchannels',badCh,'ThetaChannels',thetaCh);

%% Identify ripples 

%%%%%%%%%%% Buzcode wavelet funtion (Morlet) %%%%%%%%%%%%%%%%%%%%%%%%%%
%check where the ripple power is
whitened = WhitenSignal(lfp.data(:,126));
whitened(isnan(cleanlfp.cleantimestamps)) = nan;
msstart = 16001.000; mssend = 16011.000; % example ripple

testlfp = whitened(floor(msstart*1250):ceil(mssend*1250)); % 32



wavespec = bz_WaveSpec(wdata,'samplingRate',1250,'frange',[0.5 300],...
    'nfreqs',300,'showprogress',true,'ncyc',7);
% trade off number of frequencies and range of frequencies in order to not overwhelm memory

figure
hold all
%title(['Ch 40, ' num2str(msstart) '-' num2str(mssend) ' sec' ]);
imagesc(wavespec.timestamps,wavespec.freqs,(abs(wavespec.data))');
ax = gca;
ax.YDir = 'normal'; ax.YScale="log"; ax.YLim = [0.5 300]; ax.XLim = [min(wavespec.timestamps) max(wavespec.timestamps)];
colorbar; ax.CLim=[0 300];

%ripchan = bz_GetBestRippleChan_RW(cleanlfp); % most power in 140 to 180 Hz band. LFP has to come in as 1250 Hz
ripCh = 126;
noiseCh = 5;

% Params (defaults listed below from buzcode tutorial
% restrict periods chosen by plotting and selecting quiet period
ripdur = [30 300];% was [30 200]  %100 was default, 400 KD's suggestion
ripfreq = [140 240]; %ripfreq = [130 200]; %200 was default, 250 KD's suggestion

premaze_ms = length(premazelfp)/1250*1000;
maze_ms = length(mazelfp)/1250*1000;

wpremazelfp = WhitenSignal(lfp.data(1:lfpstartIndex-1,ripCh));
noiselfp = WhitenSignal(lfp.data(1:lfpstartIndex-1,noiseCh));
timestamps = lfp.timestamps(1:lfpstartIndex-1);
premazeRipples = bz_FindRipples(wpremazelfp,timestamps,'thresholds',[1 4],...
     'durations',ripdur,'show','on','saveMat',false,'noise',noiselfp,'EMGThresh',[]);
clear wpremazelfp noiselfp tsindex timestamps
disp('Finished finding ripples for premaze period')

wmazelfp = WhitenSignal(lfp.data(lfpstartIndex:lfpendIndex,ripCh));
noiselfp = WhitenSignal(lfp.data(lfpstartIndex:lfpendIndex,noiseCh));
timestamps = lfp.timestamps(lfpstartIndex:lfpendIndex);
mazeRipples = bz_FindRipples(wmazelfp,timestamps,'thresholds',[2 5],...
     'durations',ripdur,'show','on','saveMat',false,'noise',noiselfp,'EMGThresh',[], ...
     'restrict',[9000 max(timestamps)]); % identify wh
clear wmazelfp noiselfp tsindex timestamps
disp('Finished finding ripples for maze period')

wpostmazelfp = WhitenSignal(lfp.data(lfpendIndex+1:end,ripCh));
noiselfp = WhitenSignal(lfp.data(lfpendIndex+1:end,noiseCh));
timestamps = lfp.timestamps(lfpendIndex+1:end);
postmazeRipples = bz_FindRipples(wpostmazelfp,timestamps,'thresholds',[1 4],...
     'durations',ripdur,'show','on','saveMat',false,'noise',noiselfp,'EMGThresh',[], ...
     'restrict',[14000 18000]);
clear wpostmazelfp noiselfp tsindex timestamps
disp('Finished finding ripples for postmaze period')

ripples = vertcat(...
    [premazeRipples.timestamps(:,1) premazeRipples.peaks(:,1) premazeRipples.timestamps(:,2)], ...
    [mazeRipples.timestamps(:,1) mazeRipples.peaks(:,1) mazeRipples.timestamps(:,2)], ...
    [postmazeRipples.timestamps(:,1) postmazeRipples.peaks(:,1) postmazeRipples.timestamps(:,2)]);

SaveRippleEvents_RW(fullfile(basepath,[basename '.rip.evt']),ripples,ripCh) %ripchan is 127

ripples = ripples*1000; % to convert into ms

%% ripple triggered histogram

% 1 ) pull out all the ripple starts
% 1b) get all the spiketimes for each cell -  spiketimes in {cell}{trial}
% format
% 2 ) count number of spikes in each bin (start with 10ms) relevant to
% start of ripple
% 3) 

[spiketimes,spikeClusters,cells,cellCh,labels] = getspikeraster(phyoutputpath,...
    session.extracellular.sr,'labelstouse',[1,2,3,7,8],'showplots',true,...
'trials',trials,'plotbyfield','CleanTrial','trialsPlottype','allcells_bytrial');
% sort into which spikes belong to which cell
for c = 1:length(cells)
    spk_index = find(spikeClusters{1} == cells(c));
    c_spiketimes{c} =spiketimes{1}(spk_index);
end

keep = 1;
selectedcells = [];
for c = 1:length(cells)
    if labels(c) == 1 %|| labels(c) == 2
        selectedcells(keep,1) = keep;
       % selectedcells(2,keep) = labels(c);
        selectedcells(keep,2) = cells(c);
        keep = keep + 1;
    end
end

%%%%%%%%%%%%%%%%% fr in ripple/fr out of ripple %%%%%%%%%%%%%%%%%%%%%%%%%

%%% Sum amount of time in ms in ripple and out of ripple
inrippletime = 0;
outrippletime = 0;
lastouttime = 0;
for r = 1:length(ripples)
    intime = ripples(r,3)-ripples(r,1);
    if r == 1
        outtime(r) = ripples(r,1);
    elseif r == length(ripples)
        outtime(r) = ripples(r,1)-ripples(r-1,3);
        lastouttime = lfp.timestamps(end)*1000-ripples(r,3);
    else
        outtime(r) = ripples(r,1)-ripples(r-1,3);
    end
    inrippletime = inrippletime + intime;
    outrippletime = outrippletime + outtime(r) + lastouttime;
end

inripplefr = zeros(length(c),1); 
outripplefr = zeros(length(c),1);
inoutratio = zeros(length(c),1);

%%% number of spikes in a ripple or out
for c = 1:length(cells)
    inripple{c} = 0; outripple{c} = 0;
    lastout = 0;

    for r=1:length(ripples)
        in = c_spiketimes{c}>ripples(r,1) & c_spiketimes{c}<ripples(r,3);

        if r == 1
            out = c_spiketimes{c}<ripples(r,1);
        elseif r == length(ripples)
            out = c_spiketimes{c}<ripples(r,1) & c_spiketimes{c}>ripples(r-1,3);
            lastout = c_spiketimes{c}>ripples(r,3);
        else
            out = c_spiketimes{c}<ripples(r,1) & c_spiketimes{c}>ripples(r-1,3);
        end

        inripple{c} = inripple{c} + sum(in);
        outripple{c} = outripple{c} + sum(out) + sum(lastout);
        clear in out
    end
    inripplefr(c) = inripple{c}/inrippletime;
    outripplefr(c) = outripple{c}/outrippletime;
    inoutratio(c) = inripplefr(c)/outripplefr(c);
end

%%% plot
figure
hold all
title('inrippleSpikeRate/outrippleSpikeRate')
histogram(inoutratio,30)
xlabel('ratio')
ylabel('cell count')

%% Hilbert
fs.sr = 1250;
fs.hipass = 0.5;
fs.lopass = 500;
fs.narrowhi = 6;
fs.narrowlo = 12;
% WHITEN?:
ripCh = 126;
wmazelfp = WhitenSignal(lfp.data(lfpstartIndex:lfpendIndex,ripCh));
wlfp = WhitenSignal(lfp.data(:,ripCh));

[~,theta] = htfilter(wlfp,fs); % there'll be bad parts in here! generlaly just need to cut those
thetaphase = angle(theta);
%timestamps = lfp.timestamps(lfpstartIndex:lfpendIndex);

%%%% trial specific

plotHilbert_trials(cleanlfp,trials,fs)

for c = 1:length(cells)
    for sp = 1:length(c_spiketimes{c})
        ab = abs(lfp.timestamps*1000 - c_spiketimes{c}(sp));
        index = find(ab == min(ab));
        if length(index) > 1
            index = index(1);
        end
        c_lfpindex{c}(sp) = index;
        clear ab index
    end
    disp(['cell ' num2str(c)])
end

plotWavelet_trials(cleanlfp,trials,fs)



    c_in
   % c_index = (c_spiketimes{c}>lfp.timestamps(lfpstartIndex)*1000) .* ...
    %    (c_spiketimes{c}<lfp.timestamps(lfpendIndex)*1000); 
   % c_index = timestamps == maze_spiketimes{c};
    thphase_c{c} = thetaphase(c_spiketimes{c});
    clear c_index
end



%% divide

%%%%%%%%%%%%%%%%%%%%%% specific timing relative to peak %%%%%%%%%%%%%%%%%%

for c = 1:length(cells)
    relativeTimes{c} = [];
    for r = 1:length(ripples)
        rpeak = ripples(r,2);
        tmp = rpeak - c_spiketimes{c};
        relativeTimes{c} = vertcat(relativeTimes{c},tmp);
        clear tmp
    end
    c
end

% plot
for c = 1:7 %length(relativeTimes)
    figure
    histogram(relativeTimes{c})
    pause
end

    % rstart minus spiketimes of each cell 
%    
%         [N,edges] = histogram(relativetimes,'BinWidth',10); % 10 ms bins to start
%     end
% end
% nbins = 200;
% tbin = 0.1; % s
% rpeak = ripples(:,2);
% for c =1:length(c_spiketimes)
%     ripplecell_ccg{c} = ccg(c_spiketimes{c},rpeak,nbins,tbin);
% end
% 
% binedges = linspace(min(K),max(K),200);
% figure
% histogram(K,'BinEdges',binedges)


spks1 = c_spiketimes{2};
spks2 = c_spiketimes{3};

% spiketimes during just the maze period
for c = 1:length(cells)
    maze_spiketimes{c} = [];
    tmp = c_spiketimes{c}>lfpstartIndex & c_spiketimes{c}<lfpendIndex; 
    maze_spiketimes{c} = c_spiketimes{c}(tmp);
    clear tmp
end

%%% from CCG_tutorial.m (cell explorer)
% 
binSize = 0.01; % 1ms bin size
duration = 2; % -50ms:50ms window
[ccg,t] = CCG(maze_spiketimes,[],'binSize',binSize,'duration',duration);


figure, 
% Plotting the autocorrelogram (ACG) of the eight cell
subplot(2,1,1)
plot(t,ccg(:,3,3)), title('ASCG'), xlabel('Time (seconds)'), ylabel('Count')

% Plotting the cross correlogram (CCG) between a pair of cells
subplot(2,1,2)
plot(t,ccg(:,3,5)), title('CCG'), xlabel('Time (seconds)'), ylabel('Count')

%%% or for bar graph look
bar(t,ccg(:,3,5),'barwidth',2.5)
nplot = ceil(sqrt(size(ccg,2)));

figure
for c1 = 1:length(cells) %number of 1 rated cells %size(ccg,2)
    %hold all
    sgtitle(['Cell ' num2str(cells(c1)) ' crossCorrs'])
    for c2 = 1:30 %size(ccg,2)
        subplot(6,6,c2)
        hold on
        title(['cell ' num2str(cells(c2))])
        bar(t,ccg(:,c1,c2),'barwidth',2.5)
        hold off
    end
    pause
    clf
end





%% Explore the high frequency oscillation in cortex
%%% Using Morlet wavelet transform 
%ripplechan = 126
%%%%%%%%% Choose channel of interest and params %%%%%%%%%%%%%%%%%%%%%%%%%
corticalCh = 44; % selected by eye in neuroscope
testlfp = double(mazelfp(1000000:1500000,corticalCh)); % just pulling out the clean parts for testing
ts = lfp.timestamps(lfpstartIndex:lfpendIndex);
ts = ts(1000000:1500000);

testlfp = double(mazelfp(end-50000:end,corticalCh));

%%%%%%%% Other preliminary figures %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% power spectrum
figure
pspectrum(testlfp,1250,'FrequencyLimits',[0.5 100],'spectrogram','TimeResolution',1)

whitened = WhitenSignal(lfp.data(:,40));
whitened(isnan(cleanlfp.cleantimestamps)) = nan;
msstart = 15025.000; mssend = 15035.000;
% 818 - 828
% 4920 - 4950
% 7046 - 7056
% 7700 - 7704
% 11431 - 11441
% 13780 - 13790

testlfp = whitened(floor(msstart*1250):ceil(mssend*1250)); % 32


%%%%%%%%%%% Buzcode wavelet funtion (Morlet) %%%%%%%%%%%%%%%%%%%%%%%%%%
wavespec = bz_WaveSpec(wdata,'samplingRate',1250,'frange',[2 20],...
    'nfreqs',300,'showprogress',true,'ncyc',7);
% trade off number of frequencies and range of frequencies in order to not
% overwhelm memory

figure
hold all
hold on
%title(['Ch 126, ' num2str(msstart) '-' num2str(mssend) ' sec' ]);
imagesc(wavespec.timestamps,wavespec.freqs,(abs(wavespec.data))');
ax = gca;
ax.YDir = 'normal';  ax.YLim = [2 20]; ax.XLim = [min(wavespec.timestamps) max(wavespec.timestamps)];
colorbar; ax.CLim=[0 300];ax.YScale="log"; ax.xLabel = 'Time (s) from trial start'; 
ax.yLabel = 'Log Frequency';

%%%%%%%%%%%%%%%%%% Get low res, wide freq figure %%%%%%%%%%%%%%%%%%%%%%%

%%% For maze, in chunks
winsize = 300; % 5 minutes
For w = 1
wavespec = bz_WaveSpec(testlfp,'samplingRate',1250,'frange',[1 128],...
    'nfreqs',128,'showprogress',true,'ncyc',7);

%%% For post, in chunks


%%% For pre, in chunks


%%%%%%%%%%%%%%%%%% Get high res, theta only figure %%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%% Get high res, high freq only figure %%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%% Output results %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Also post it notes

%%%%%%%%%%%%%%%%%% Narrowing in on specific band of interest %%%%%%%%%%

testlfp = cleanlfp.data(floor(6156833/1000*1250):ceil(6157833/1000*1250),39:49);

params.Fs = 1250;
params.fpass = [20 40];
[S,t,f] = mtspecgramc(testlfp,[1 0.5],params);

figure
imagesc(t,f,S)

%%% Looking at individual units with timestamps in a ripple period %%%%%%%%
cstart = 3729700; % closest points are 3729000 and 3730000? 
cstop = 3729710;
inrippleSpikes = spiketimes{1}>cstart & spiketimes{1}<cstop;
idx = find(inrippleSpikes == 1);
inrippleCells = unique(spikeClusters{1}(idx));
allinrippleCells = spikeClusters{1}(idx);

keep = 1;
for c = 1:length(allinrippleCells)
    tmp = find(cells == allinrippleCells(c));
    if isempty(tmp) % labeled zero
        continue
    else
        cellIdx(keep) = tmp;
        keep = keep + 1;
    end
end

inrippleLabels = labels(cellIdx); % 3 7 1 for the first test case

%% Adapting Utku HighVoltageSignals code %%%%%%%%%%%%%%%%%%

% get the spectrogram first
params.Fs = 1250;
params.fpass = [1 40];
[S,t,f] = mtspecgramc(testlfp,[2,1],params);

%%%%% get timestamps for high voltage signals 

corticalCh = 39:53;
hpcCh = 119:127;

premaze_ms = length(premazelfp)/1250*1000;
maze_ms = length(mazelfp)/1250*1000;

wpremazelfp = WhitenSignal(lfp.data(1:lfpstartIndex-1,corticalCh));
remove = isnan(cleanlfp.cleantimestamps(1:lfpstartIndex-1));
wpremazelfp(remove,:) = nan;
hvs_premaze = getHighVoltageSignals(wpremazelfp,corticalCh,sessionInfo.lfpSampleRate);
hvs_premaze = hvs_premaze*1000;
clear wpremazelfp remove
disp('Finished finding HVS from premaze period')

wmazelfp = WhitenSignal(lfp.data(lfpstartIndex:lfpendIndex,corticalCh));
remove = isnan(cleanlfp.cleantimestamps(lfpstartIndex:lfpendIndex));
wmazelfp(remove,:) = nan;
hvs_maze = getHighVoltageSignals(wmazelfp,corticalCh,sessionInfo.lfpSampleRate);
hvs_maze = premaze_ms + hvs_maze*1000;
clear wmazelfp remove
disp('Finished finding HVS from maze period')

wpostmazelfp = WhitenSignal(lfp.data(lfpendIndex+1:end,corticalCh));
remove = isnan(cleanlfp.cleantimestamps(lfpendIndex+1:end));
wpostmazelfp(remove,:) = nan;
hvs_postmaze = getHighVoltageSignals(wpostmazelfp,corticalCh,sessionInfo.lfpSampleRate);
hvs_postmaze = premaze_ms + maze_ms + hvs_postmaze*1000;
clear wpostmazelfp remove
disp('Finished finding HVS from postmaze period')

hvsEvents = vertcat(hvs_premaze, hvs_maze, hvs_postmaze);

savefilename = [basepath '/' basename '.hvs.evt'];
channelIDs = '39:53';
SaveHVSEvents_RW(savefilename,hvsEvents,channelIDs)

%% Determine VTE like behavior
% try the Papale code again, or whiteboard your own definition usin the
% behavior rotational data


%%%%%%%%%%%%%%%%%% White board what makes sense %%%%%%%%%%%%%%%%%%%%%%%%%
% 1) get position data, determine if they're in the choice box
% 2) get rotation data for that time period
% 3) get trial direction; set the final rotation orientation
% 4) OBSERVE what the rotation looks like before then; see if there's a
% pattern
%   a) hypothesis is that if there is a rotation in the opposite direction pre
%   final direction, it's a VTE. however, the rat was stationary a lot, so
%   it might be harder to define for this animal

boundaries = sections_RW(optitrack.position.interpolatedx,optitrack.position.interpolatedy);
newtrials = identifyVTE(trials,boundaries,'checkplot',true);

%% Sorting spikes into spatial bins
% Adam code

%%%%%%%%%%%%%% Run it as if it were HPC %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tuningcurves2D(behaviortimestamps,phyoutputpath,'trials',trials,'labelstouse',1:2,'showplots',true); %false);

[FRmaps,xbins,ybins,empties] = place_cells03_RW(behaviortimestamps,optitrack,phyoutputpath,'labelstouse',1:2,'placecellplots',true); %false);

%%%%%%% ID if any cells are HPC? Just try running those to compare %%%%%%%
%% Do some cells fire more on L vs R trials? 
% adapt the code you've already almost got 

%%%%%%%%%%% Adapt the trial specific code %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nleft = length(find(strcmp(trials.trialDirection,'Left')==1));
nright = length(find(strcmp(trials.trialDirection,'Right')==1));

% plots and saves plot
getspikeraster(phyoutputpath,session.extracellular.sr,...
    'trials',trials,'trialsPlottype','condition_withincells',...
    'plotbyfield','trialDirection')


%% F31 example figure

% 3 by 2
% (3,2,1) = all behavioral trials
% (3,2,2) = place cell fig
% (3,2,3:4) = bandpassed theta + gamma
% (3,2,5:6) = raster for all cells one trial (same trial as th/ga)

tri = 6;

%%% 
subplot(3,2,2)
hold all
for t= 1:size(trials,1)
    plot(-trials.ytrialdata{t},trials.xtrialdata{t});
end
title('All trials for sample day')


subplot(3,2,5:6)
tstart_s = trials.startTrial_file(tri)/1000; % into seconds
tend_s = trials.endTrial_file(tri)/1000;
hold on
plotRaster_resClu(data,[tstart_s tend_s]);

%%%%%%%%%%% Buzcode wavelet funtion (Morlet) %%%%%%%%%%%%%%%%%%%%%%%%%%
%check where the ripple power is

dif = lfptimestamps.fromlfpfilestart - trials.startTrial_file(6);
tstart_idx = find(abs(dif) == min(abs(dif)));
dif = lfptimestamps.fromlfpfilestart - trials.endTrial_file(6);
tend_idx = find(abs(dif) == min(abs(dif)));

whitened = WhitenSignal(lfp.data(tstart_idx:tend_idx,126));
whitened(isnan(cleanlfp.cleantimestamps)) = nan;

wavespec = bz_WaveSpec(whitened,'samplingRate',1250,'frange',[5 12],...
    'nfreqs',300,'showprogress',true,'ncyc',7);
% trade off number of frequencies and range of frequencies in order to not overwhelm memory

subplot(3,2,5:6)
hold all
%title(['Ch 40, ' num2str(msstart) '-' num2str(mssend) ' sec' ]);
imagesc(wavespec.timestamps,wavespec.freqs,(abs(wavespec.data))');
ax = gca;
ax.YDir = 'normal'; ax.YScale="log"; ax.YLim = [5 12]; ax.XLim = [min(wavespec.timestamps) max(wavespec.timestamps)];
ax.CLim=[0 300]; colorbar; 
ylabel('Frequency (Hz)');


%%%%%%%%%%%%%%%%%%%% UK code %%%%%%%%%%%%%%%%%%%%%%%%%%%

sf = neuro.spike.SpikeFactory.instance();
sa = sf.getPhyOutputFolder;

optifolder = '/data/ExperimentsRSW/CircularMaze/20230614/Optitrack';
op = optiTrack.OptiLoader.instance(optifolder);
op1 = op.loadFile; %will open the folder for selection. do op1, op2, op3 .. for all files
op2 = op.loadFile;
op3 = op.loadFile;
op4 = op.loadFile;
opti = optiTrack.OptiFileCombined(op1,op2,op3,op4);
positionData = opti.getMergedPositionData;

%to visualize data:
%plot(positionData.data.X,positionData.data.Z,'LineWidth',2)










