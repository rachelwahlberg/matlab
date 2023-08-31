function [behaviortimestamps,lfptimestamps] = getRelativeTimestamps(varargin)
% Rachel Wahlberg 7/22
% Aligning behavior files and lfp files

p = inputParser;

addRequired(p,'basepath',@ischar)
addRequired(p,'behaviorpath',@ischar)
addParameter(p,'baseLFPfolders',[],@iscell); %check if cell
addParameter(p,'TimeIntervalCombinedpath',[],@ischar) %path to time interval combined (utku merge file)
addParameter(p,'rawlfpsamprate',30000,@isnumeric)
addParameter(p,'showplot',false,@islogical)
addParameter(p,'savelfpfilename',[],@ischar)
addParameter(p,'savebehaviorfilename',[],@ischar)

parse(p,varargin{:})
basepath = p.Results.basepath;
basename = bz_BasenameFromBasepath(p.Results.basepath);
behaviorpath = p.Results.behaviorpath;
baseLFPfolders = p.Results.baseLFPfolders;
TimeIntervalCombinedpath = p.Results.TimeIntervalCombinedpath;
rawlfpsamprate = p.results.rawlfpsamprate;
savelfpfilename = p.Results.savelfpfilename;
savebehaviorfilename = p.Results.savebehaviorfilename;
showplot = p.Results.showplot;

%% Set basepaths/directories
%basepath = '/home/wahlberg/Exp_Data/M1_Nov2021/20211123/';
%basename = 'merged_M1_20211123_raw';
%behdir = dir(fullfile(basepath,'Behavioral/files_to_use','*.csv'));

%% Get starttimes from lfp files in datetime format
% RAW lfp timestamps so sampling rate is 30000 not 1250! 
if isfile([basename '.lfptimestamps.mat']) % '_1250Hz_lfptimestamps.mat'
    load([basename '.lfptimestamps.mat'])
    disp(['Loading ' basename '.lfptimestamps.mat'])
else

    clear fID dataArray
    formatSpec = '%q%[^\n\r]';
    % '/home/wahlberg/Exp_Data/M1_Nov2021/20211123/2021-11-23_15-34-56/RecordNode108',...
    if ~isempty(baseLFPfolders)
        i = 1;
        for folder = 1:length(baseLFPfolders)
            d = dir(fullfile(baseLFPfolders{folder}, 'settings*'));

            fID = fopen(fullfile(baseLFPfolders{folder},d(1).name),'r');
            dataArray = textscan(fID,formatSpec,'Delimiter', ',', 'ReturnOnError', false);

            tmp = split(dataArray{1}{6});
            tmp = split(tmp{4},'<');
            lfpfile.starttime{i} = tmp{1};
            i = i+1;

            if length(d) > 1 % if two experiments
                clear fID dataArray
                fID = fopen(fullfile(baseLFPfolders{folder},d(2).name),'r');
                dataArray = textscan(fID, formatSpec, 'Delimiter', ',',  'ReturnOnError', false);

                tmp = split(dataArray{1}{6});
                tmp = split(tmp{4},'<');
                lfpfile.starttime{i} = tmp{1};
                i = i+1;
            end
            clear fID d dataArray
        end
    end

    % shouldn't the mergedrelative time just be equal to the original
    % timepoints?
    if TimeIntervalCombinedpath
        % TimeIntervalCombined is file from utku's merging files code
        formatSpec = '%q%q%q%[^\n\r]';
        fID = fopen(TimeIntervalCombinedpath,'r');
        dataArray = textscan(fID, formatSpec, 'Delimiter', ',',  'ReturnOnError', false);

        n = 1;
         lfpfile.mergedrelativetime = timetable(seconds(0));
        for i = 4:length(dataArray{2}) % 4 because not using some of the earlier empty files
            lfpfile.nsmps(n) = str2double(dataArray{2}{i});
            lfpfile.starttime{n} = datetime(dataArray{1}{i},'inputFormat','yyyy-MM-dd HH:mm:ss.SSS','Format','yyyy-MM-dd HH:mm:ss.SSS');
            s = linspace(0,lfpfile.nsmps(n)/rawlfpsamprate,lfpfile.nsmps(n));
            s = seconds(s);

            lfpfile.realtime{n} = s + lfpfile.starttime{n}; % actual time of day
            lfpfile.relativetime{n} = lfpfile.realtime{n} - lfpfile.starttime{1}; % time relative to start of first lfpfile

            reltime = timetable((lfpfile.relativetime{n})');
            lfpfile.mergedrelativetime = [lfpfile.mergedrelativetime;reltime];
            
            n = n+1;
        end

        nlfp = length(lfpfile.nsmps); % total number of lfp files

        % put into ms from start of the lfp file
        lfpfile.mergedrelativetime = seconds(lfpfile.mergedrelativetime.Time(:)*1000);
        lfpfile.mergedrelativetime(1) = []; % to get rid of the original zero allocation
    end

    if isempty(baseLFPfolders) && isempty(TimeIntervalCombinedpath)
        disp('Specify where to pull LFP filetimes from')
    end

end

%% Get starttimes from behavior in datetime format
% Have folder with all the csv files from optitrack

if isfile([basename '.behaviortimestamps.mat']) %'_60hz_behaviortimestamps.mat'
    load([basename '.behaviortimestamps.mat'])
    disp(['Loading ' basename '.behaviortimestamps.mat'])
else

    behdir = dir(fullfile(behaviorpath,'*.csv'));

    % should have folder "behaviorfiles" in basepath
    formatSpec = '%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%[^\n\r]';
    header_length = 7;

    behaviortimestamps.mergedrelativetime = timetable(seconds(0));
    for f = 1:length(behdir)
        fID = fopen(fullfile(behaviorpath, behdir(f).name),'r');
        dataArray = textscan(fID, formatSpec, 'Delimiter', ',',  'ReturnOnError', false);

        behaviortimestamps.starttime{f} = datetime(dataArray{12}{1},'inputFormat','yyyy-MM-dd hh.mm.ss.SSS a','Format','yyyy-MM-dd HH:mm:ss.SSS');
        s = (seconds(str2double(dataArray{2}(header_length:end))))';
        behaviortimestamps.realtime{f} = s + behaviortimestamps.starttime{f};
        behaviortimestamps.relativetime{f} = behaviortimestamps.realtime{f}-lfpfile.starttime{1}; %relative to start of first LFP file  
        behaviortimestamps.nsmps(f) = str2double(dataArray{16}{1}); % total frames

        reltime = timetable((behaviortimestamps.relativetime{f})');
        behaviortimestamps.mergedrelativetime = [behaviortimestamps.mergedrelativetime;reltime];

        behaviortimestamps.framerate{f} = dataArray{8}{1};
        clear dataArray
    end

    % put into MILLISECONDS relative to lfp start, get into double format
    behaviortimestamps.mergedrelativetime(1,:) = []; % to get rid of the original zero allocation
    behaviortimestamps.mergedrelativetime = seconds(behaviortimestamps.mergedrelativetime(:)*1000);

    nbeh = length(behaviortimestamps.starttime); % total number of behavior files
end

lfptimestamps.notes.alignment = 'time 0 is the starttime of the first LFP file';
lfptimestamps.notes.samplingrate = lfpsamplingrate;

behaviortimestamps.notes.alignment = 'time 0 is the starttime of the first LFP file';
behaviortimestamps.notes.samplingrate = behaviorsamplingrate;

% include notes about conversion, original sampling rates, etc
%% Check that it looks right! Also you can visually check if any files seem to be missing

if showplot == 1
    figure
    hold all
    scatter(behaviortimestamps.mergedrelativetime.Time,1:length(behaviortimestamps.mergedrelativetime.Time))
    scatter(lfpfile.mergedrelativetime.Time,1:length(lfpfile.mergedrelativetime.Time),'r')
end

%% Save

if savelfpfilename
    save(savelfpfilename,'lfpfile')
end

if savebehaviorfilename
    save(savebehaviorfilename,'behaviortimestamps')
end

%% Other checks

% both of these should be about equal and should make sense for the length of time spent recording
% these are how you'll align the times
%behavior.reltime_merged.Time(end)
%lfpfile.reltime_merged.Time(end)
%% Load data to add in timeframes
%
% basepath = '/home/wahlberg/Exp_Data/M1_Nov2021/20211123/merged_M1_20211123_raw';
% optitrack = optitrack2buzcode('basepath',basepath,'basename','merged_M1_20211123_raw',...
%     'filenames',{'20211123_03.33.38PM.csv','20211123_03.44.00PM.csv', ...
%     '20211123_03.51.19PM.csv','20211123_04.00.42PM.csv','20211123_04.20.37PM.csv','20211123_04.41.49PM.csv'},saveMat=false);
%
%
% lt = lfp_timestamps.reltime_merged.Time;
% bt = behavior_timestamps.reltime_merged.Time;
% behavior_timestamps.timestamps = seconds(bt - lt(1));
% lfp_timestamps.timestamps = seconds(lt-lt(1));
%
% % these are the already merged datapoints?
%
% lfp = bz_GetLFP('all','basepath',fullfile(basepath,basename)); % or use m =memmapfile(merged_M1_20211123_raw.lfp','Format','int16') and then do the channels/division by SR manually
% optitrack = load(fullfile(basepath,basename,[basename '.optitrack.behavior.mat']));
%
%
%
% interp_vec = linspace(1,length(lfp.timestamps)




