%% AnalysisPipeline
% Mostly Cell Explorer/Buzcode, with additions from myself/
% Diba lab/Brendon Watson

%% 0a) merge files
% Run only pre-sorting! 

usedate = 112421;
basefolders = getbasefolders(usedate); % adapt to include the variables below to not have to rerun each time
%mergedfilename = '/home/wahlberg/Analysis/M1_Nov2021/20211124/merged_M1_20211124';
mergedfilename = '/data/20211124/merged_M1_20211124_test';
OpenEphysPreprocess_RW(basefolders,mergedfilename);

%%% 0b) spyking circus in merged folder - need .params, .prb, .dat, dead.txt (if bad periods)
%%% 0c) phy in merged folder

%%%%%%%%%%%%%%%%%% 1) Define paths, session info, timestamps %%%%%%%%%%%%%%
%% 1a) Define the basepath of the dataset to run. 
% The dataset should at minimum consist of the raw data and spike sorted data.
basepath = '/data/20211124/merged_M1_20211124';
basename = 'merged_M1_20211124';
behaviorpath = fullfile(basepath, 'behaviorfiles/');
%manualbehaviorpath = '/home/wahlberg/Exp_Data/M1_Nov2021/20211123/Behavioral/221123_manualscoring.csv';
phyoutputpath =  fullfile(basepath,'merged_M1_20211124/merged_M1_20211124crs.GUI');
% NEED THIS OR
% baseLFPfolders={'/home/wahlberg/Analysis/M1_Nov2021/20211123/2021-11-23_15-51-17/RecordNode108',...
%     '/home/wahlberg/Analysis/M1_Nov2021/20211123/2021-11-23_16-00-39/RecordNode108',...
%     '/home/wahlberg/Analysis/M1_Nov2021/20211123/2021-11-23_16-20-34/RecordNode108',...
%     '/home/wahlbesession = sessionTemplate(pwd,'showGUI',true); % CELL EXPLORER
%rg/Analysis/M1_Nov2021/20211123/2021-11-23_16-41-47/RecordNode108'};
% THIS
TimeIntervalCombinedpath = fullfile(basepath, [basename '_1250Hz.TimeIntervalCombined.csv']);

cd(basepath)
%% 1b) CE Generate session metadata struct using the template function and display the meta data in a gui

session = sessionTemplate(pwd,'showGUI',false); % CELL EXPLORER
validateSessionStruct(session); % validate the required and optional fields

sessionInfo = bz_getSessionInfo(basepath); % BUZCODE see < https://github.com/buzsakilab/buzcode/wiki/Data-Formatting-Standards > . 

% You can also inspect the session struct by calling the GUI directly:
% session = gui_session(session);

%%%%%%%%%%%%%%%%%% 2) Get LFP and behavior files %%%%%%%%%%%%%%%%%%%%%%%%%%
%% 2a) BZ Get the lfp data from the dat file
% lfp is a struct
if isfile(fullfile(basepath, [basename '.lfp']))
    disp('Loading existing LFP file')
    lfp = bz_GetLFP('all','basepath',basepath,'basename',basename);
    %load(fullfile(basepath, [basename '_lfp.mat']));
else
    disp('LFP file does not yet exist, creating one from memmap')
    % unless 
    %m =memmapfile('merged_M1_20211124.lfp','Format','int16'); % and then do the channels/division by SR manually
    %lfp2 = (reshape(m.Data,session.extracellular.nChannels,[]))';

    %bz_LFPfromDat(basepath) % DON'T USE, doesn't work correctly!
    lfp = bz_GetLFP('all','basepath',basepath,'basename',basename); 
    % or use m =memmapfile(merged_M1_20211123_raw.lfp','Format','int16') and then do the channels/division by SR manually
end
% lfp.timestamps is in seconds
%lfp.timestamps = lfp.timestamps*1000;
%% 2b) BZ/RW Import behavioral tracking from optitrack
% first we import the raw behavioral tracking into a Matlab struct:
% the optitrack output contains timestamps (from optitrack) and position data 
%    optitrack.timestamps      timestamps from optitrack
%    optitrack.position.x      x-position
%    optitrack.position.y      y-position
%    optitrack.position.z      z-position
%    optitrack.speed           speed
%    optitrack.sr              samplingrate

% then add interpolated data points and rotate if necessary

if isfile(fullfile(basepath, [basename '.dat.optitrack.behavior.mat']))
    disp('Loading existing optitrack file')
    load(fullfile(basepath, [basename '.dat.optitrack.behavior.mat']));
else
    optitrack = getbehavior(basepath,basename,behaviorpath);
end

%% 2c) RW Get timestamps for artifacts/remove bad channels

% Due to memory issues! the buzcode loaded LFP is too big for longer recordings

[cleanlfp,removetimestamps] = glitchdetector(lfp,session.extracellular.srLfp,...
    'percentiles',[0.9 99.93]); %,'saveMat',true,'eventfilename',fullfile(basepath,[basename '.bad.evt']));

%for s= 1:length(removetimestamps.sampformat)

%testlfp = 
%% 2d) Get timestamps of LFP and behavior files relative to start of LFP
% periods with "removetimestamps" will be nans if included

if isfile(fullfile(basepath, [basename '.lfptimestamps.mat'])) && isfile(fullfile(basepath, [basename '.behaviortimestamps.mat']))
    load(fullfile(basepath, [basename '.lfptimestamps.mat']));
    load(fullfile(basepath, [basename '.behaviortimestamps.mat']));
else
    [lfptimestamps,behaviortimestamps] = getTimestamps(basepath,behaviorpath,...
        'TimeIntervalCombinedpath',TimeIntervalCombinedpath,'removetimestamps',removetimestamps);
end
%startofbehavior = behtimestamps.reltime_merged.Time(1); % startofbehavior relative to start of LFP, in ms. ADD this value

%%%%%%%%%%%%%%%%%% 3) Divide into trials %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 3a) Divide into trials
if isfile([basename '.trials.mat'])
    load([basename '.trials.mat'])
else
    newtrialfilename = 'basename.trials.mat';
trials = gettrials(optitrack,behaviortimestamps,'alternationYmaze','savetrials',true,'plottrials',true);
% in clock ms from start of LFP file
end

%%%%%%%%%%%%%%%%%% 4) Get spike rasters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Get spike rasters for left vs right trials

nleft = length(find(strcmp(trials.trialDirection,'Left')==1));
nright = length(find(strcmp(trials.trialDirection,'Right')==1));

% plots and saves plot
getspikeraster(phyoutputpath,session.extracellular.sr,...
    'trials',trials,'plotbyfield','trialDirection')

%% Get place fields
% clustering_path = '/home/wahlberg/Exp_Data/M1_Nov2021/20211123/merged_M1_20211123/merged_M1_20211123/merged_M1_20211123crs.GUI';
% Phy2Neurosuite(basepath,clustering_path)
% spikes = bz_GetSpikes('basepath',basepath,'saveMat',false);

%%% 2D firing rate map
%placecellplots = 1 will pause between each cell
[FRmaps,xbins,ybins,empties] = place_cells03_RW(behaviortimestamps,optitrack,phyoutputpath,'labelstouse',1:3,'placecellplots',true); %false);

%%%%%%%%%%%%%%%%%% 5) Get events of interest %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 4a) BZ Get SWR
% doesn't work with super noisy data, so using "cleanlfp" - 
% RW version removes nans from cleanlfp in order to calculate best channel
ripchan = bz_GetBestRippleChan_RW(cleanlfp); % most power in 140 to 180 Hz band. LFP has to come in as 1250 Hz

ripchan = 126;
% Params (defaults listed below from buzcode tutorial
ripthresh = [1 4]; % was [3 5] %min std, max std [2 5] default
ripdur = [30 400];% was [30 200]  %100 was default, 400 KD's suggestion
ripfreq = [130 200]; %ripfreq = [130 200]; %200 was default, 250 KD's suggestion

riplfp = cleanlfp.data(:,ripchan);
riplfp(isnan(riplfp)) = [];

noiselfp = cleanlfp.data(:,5);
noiselfp(isnan(noiselfp)) = [];

timestamps = cleanlfp.cleantimestamps; %
timestamps(isnan(cleanlfp.cleantimestamps)) = [];

% doing a super crappy job (missing many ripples and identifying many false ripples)
ripples = bz_FindRipples(riplfp,timestamps,'thresholds',ripthresh,...
     'durations',ripdur,'show','on','saveMat',false,'noise',noiselfp,'EMGThresh',[]);

% doing a super crappy job (missing many ripples and identifying many false ripples)
% try same code again with data with stronger ripple power in the first
% place!
ripples = ce_FindRipples(riplfp,timestamps,'thresholds',ripthresh,...
    'durations',ripdur,'show','on','saveMat',false,'noise',noiselfp,'EMGThresh',[]);
SaveRippleEvents_RW(fullfile(basepath,[basename '.rip.evt']),ripples,ripchan) %ripchan is 127

% % Calculate ripple stats (not currently working)
% % then sort by duration vs amplitude of ripple
% [maps,data,stats] = bz_RippleStats(double(riplfp),timestamps,ripples);
% [~,dursort]=sort(data.duration,1,'descend');
% [~,ampsort]=sort(data.peakAmplitude,1,'descend');

%% Get sleep states
% Aug 2022 not working great - but because there's no sleep in the test
% recording?
badCh = session.channelTags.Bad.channels;
thetaCh = session.channelTags.Theta.channels; % check this
SleepState = SleepScoreMaster(basepath,'ignoretime',removetimestamps.sampformat,'rejectchannels',badCh,'ThetaChannels',thetaCh);

%% CE 3.1.1 Run the cell metrics pipeline 'ProcessCellMetrics' using the session struct as input
cell_metrics = ProcessCellMetrics('session', session,'showGUI',false);

%OR %%
cell_metrics = load("merged_M1_20211123_raw.dat.cell_metrics.cellinfo.mat");

%% CE 3.1.2 Visualize the cell metrics in CellExplorer
cell_metrics = CellExplorer('metrics',cell_metrics); 

%% CE 3.2 Open several session from basepaths
basepaths = {'/your/data/path/basename_1/','/your/data/path/basename_2/'};
cell_metrics = loadCellMetricsBatch('basepaths',basepaths);
cell_metrics = CellExplorer('metrics',cell_metrics);

%% BZ 3. Intro spectral processing: filter in ripple frequency and compare to output of bz_FindRipples
% 
% % Filter channel in ripple freq range
%     ripfiltlfp = bz_Filter(lfp.data,'passband',ripfreq,'filter','fir1'); % lots of options in bz_Filter, read through it
% 
%     %% BZ 4.   Plots
%     %% BZ 4i.  Look: Plot some ripples
%     figure()
%     x=maps.ripples(dursort,:);
%     for i=1:100
%         subplot(10,10,i)
%         plot(x(i,:))
%         set(gca,'XTick',[]);
%         set(gca,'YTick',[]);
%         axis off
%     end

%% Try out CMBHOME (Hasselmo code)



spk(1,1) = CMBHOME.Spike('ts', my_spike_times_t1_c1);
spk(1,2) = CMBHOME.Spike('ts', my_spike_times_t1_c2);
spk(12,3) = CMBHOME.Spike('ts', my_spike_times_t12_c3);

lfp(1) = CMBHOME.LFP(t1_sig,t1_ts,t1_fs);
lfp(12) = CMBHOME.LFP(t12_sig, t12_ts, t12_fs);

root = CMBHOME.Session('b_x', x, 'b_y', y, ...
                       'b_ts', ts, 'b_headdir', hd, ...
                       'fs_video', fs_video, ...
                       'spike', spk, ...
                       'b_lfp', lfp, ...
                       );

labelstouse = 1:2;
spikes = getspikes(phyoutputpath,'labelstouse',labelstouse); %,'lfptimestamps',lfptimestamps);
% get the spikes into logical format
cmb_spk(1,1) = CMBHOME.Spike('ts',spikes.logical);
cmb_lfp(1) = CMBHOME.LFP(lfp.data,lfptimestamps.fromlfpfilestart,lfp.samplingRate);



%% EEG states

testlfp2 = testlfp(1:1000000,1);
CheckEegStates_RW(basepath,basename,'lfp',testlfp2,'redo_flag',true);
















