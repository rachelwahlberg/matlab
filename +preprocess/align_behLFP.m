function [] = align_behLFP(varargin)
% Rachel Wahlberg 7/22
% Aligning behavior files and lfp files

p = inputParser;

addRequired(p,'basepath',@isstr)
addRequired(p,'behaviorpath',@isstr)
addParameter(p,'saveMat',false,@islogical);
addParameter(p,'baseLFPfolders',[],@iscell); %check if cell
addParameter(p,'TimeIntervalCombinedpath',[],@isstr) %path to time interval combined (utku merge file)

parse(p,varargin{:})
basepath = p.Results.basepath;
basename = bz_BasenameFromBasepath(p.Results.basepath);
behaviorpath = p.Results.behaviorpath;
baseLFPfolders = p.Results.baseLFPfolders;
TimeIntervalCombinedpath = p.Results.TimeIntervalCombinedpath;
saveMat = p.Results.basepath;


%% Set basepaths/directories
%basepath = '/home/wahlberg/Exp_Data/M1_Nov2021/20211123/';
%basename = 'merged_M1_20211123_raw';
%behdir = dir(fullfile(basepath,'Behavioral/files_to_use','*.csv'));

%% Get starttimes from lfp files in datetime format
clear fID dataArray
% formatSpec = '%q%[^\n\r]';
%
% '/home/wahlberg/Exp_Data/M1_Nov2021/20211123/2021-11-23_15-34-56/RecordNode108',...
if baseLFPfolders
    i = 1;
    for folder = 1:length(basefolders)
        d = dir(fullfile(basefolders{folder}, 'settings*'));

        fID = fopen(fullfile(basefolders{folder},d(1).name),'r');
        dataArray = textscan(fID, formatSpec, 'Delimiter', ',',  'ReturnOnError', false);

        tmp = split(dataArray{1}{6});
        tmp = split(tmp{4},'<');
        lfpfile.starttime{i} = tmp{1};
        i = i+1;

        if length(d) > 1 % if two experiments
            clear fID dataArray
            fID = fopen(fullfile(basefolders{folder},d(2).name),'r');
            dataArray = textscan(fID, formatSpec, 'Delimiter', ',',  'ReturnOnError', false);

            tmp = split(dataArray{1}{6});
            tmp = split(tmp{4},'<');
            lfpfile.starttime{i} = tmp{1};
            i = i+1;
        end
        clear fID d dataArray
    end
end

if TimeIntervalCombinedpath

    % Alternatively you could've just gotten it from time interval combined
    % file too (from Utku's merging files code)
    formatSpec = '%q%q%q%[^\n\r]';
    fID = fopen(fullfile(basepath,basename,[basename '_1250Hz.TimeIntervalCombined.csv']),'r');
    dataArray = textscan(fID, formatSpec, 'Delimiter', ',',  'ReturnOnError', false);

    n = 1;
    m = 0;
    %lfpfile.reltime_merged = 0; %on purpose zero instead of [] because the merged file is one sample greater than the sum of the parts
    lfpfile.reltime_merged = timetable(seconds(0));
    for i = 4:length(dataArray{2}) % 4 because not using some of the earlier empty files
        lfpfile.nsmps(n) = str2double(dataArray{2}{i}); % at 1250 Hz
        lfpfile.starttime{n} = datetime(dataArray{1}{i},'inputFormat','yyyy-MM-dd HH:mm:ss.SSS','Format','yyyy-MM-dd HH:mm:ss.SSS');
        s = linspace(0,lfpfile.nsmps(n)/1250,lfpfile.nsmps(n));
        s = seconds(s);
        % lfpfile.reltime_merged = cat(2,lfpfile.reltime_merged,s);
        lfpfile.realtime{n} = s + lfpfile.starttime{n}; % actual time of day
        lfpfile.relativetime{n} = lfpfile.realtime{n} - lfpfile.starttime{1}; % time relative to start of first lfpfile

        reltime = timetable((lfpfile.relativetime{n})');
        lfpfile.reltime_merged = [lfpfile.reltime_merged;reltime];

        % lfpfile.reltime_merged(m+1:m+lfpfile.nsmps(n)) = lfpfile.relativetime{n}; % relative time

        m = m+lfpfile.nsmps(n);
        n = n+1;
    end

    nlfp = length(lfpfile.nsmps); % total number of lfp files

end

if isempty(baseLFPfolders) && isempty(TimeIntervalCombinedpath)
    disp('Specify where to pull LFP filetimes from')
end

%% Get starttimes from behavior in datetime format
% Have folder with all the csv files from optitrack
behdir = dir(fullfile(behaviorpath,'*.csv'));

% should have folder "behaviorfiles" in basepath
formatSpec = '%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%[^\n\r]';
header_length = 7;

m = 0;
behavior.reltime_merged = timetable(seconds(0));
for f = 1:length(behdir)
    fID = fopen(fullfile(behaviorpath, behdir(f).name),'r');
    dataArray = textscan(fID, formatSpec, 'Delimiter', ',',  'ReturnOnError', false);

    behavior.starttime{f} = datetime(dataArray{12}{1},'inputFormat','yyyy-MM-dd hh.mm.ss.SSS a','Format','yyyy-MM-dd HH:mm:ss.SSS');
    s = (seconds(str2double(dataArray{2}(header_length:end))))';
    behavior.realtime{f} = s + behavior.starttime{f};
    behavior.relativetime{f} = behavior.realtime{f}-lfpfile.starttime{1}; %relative to start of first LFP file
    behavior.nsmps(f) = str2double(dataArray{16}{1}); % total frames

    %behavior.reltime_merged = cat(2,behavior.reltime_merged,s);
    %behavior.reltime_merged(m+1:m+behavior.nsmps(f)) = behavior.relativetime{f}; % relative time

    reltime = timetable((behavior.relativetime{f})');
    behavior.reltime_merged = [behavior.reltime_merged;reltime];

    behavior.framerate{f} = dataArray{8}{1};
    m = m+behavior.nsmps(f);
    clear dataArray
end
behavior.reltime_merged(1,:) = []; % to get rid of the original zero allocation

nbeh = length(behavior.starttime); % total number of behavior files

%% Check that it looks right! Also you can visually check if any files seem to be missing

if showplot == 1
figure
hold all
scatter(behavior.reltime_merged.Time,1:length(behavior.reltime_merged.Time))
scatter(lfpfile.reltime_merged.Time,1:length(lfpfile.reltime_merged.Time),'r')

end

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







































