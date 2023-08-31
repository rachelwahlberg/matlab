function [lfptimestamps,behaviortimestamps] = getTimestamps(varargin)
% Rachel Wahlberg 7/22
% Aligning behavior files and lfp files

p = inputParser;

addRequired(p,'basepath',@ischar)
addRequired(p,'behaviorpath',@ischar)
addParameter(p,'baseLFPfolders',[],@iscell); %check if cell
addParameter(p,'TimeIntervalCombinedpath',[],@ischar) %path to time interval combined (utku merge file)
addParameter(p,'rawlfpsamprate',30000,@isnumeric)
addParameter(p,'removetimestamps',[],@isstruct) % with field .sampformat [start end], and .sampfreq for samples
addParameter(p,'showplot',false,@islogical)
addParameter(p,'savelfpfilename',[],@ischar)
addParameter(p,'savebehaviorfilename',[],@ischar)

parse(p,varargin{:})
basepath = p.Results.basepath;
basename = bz_BasenameFromBasepath(p.Results.basepath);
behaviorpath = p.Results.behaviorpath;
baseLFPfolders = p.Results.baseLFPfolders;
TimeIntervalCombinedpath = p.Results.TimeIntervalCombinedpath;
rawlfpsamprate = p.Results.rawlfpsamprate;
removetimestamps = p.Results.removetimestamps;
savelfpfilename = p.Results.savelfpfilename;
savebehaviorfilename = p.Results.savebehaviorfilename;
showplot = p.Results.showplot;

%% Set basepaths/directories
%basepath = '/home/wahlberg/Exp_Data/M1_Nov2021/20211123/';
%basename = 'merged_M1_20211123_raw';
%behdir = dir(fullfile(basepath,'Behavioral/files_to_use','*.csv'));

%% Get starttimes from lfp files in datetime format

%realtime is in ms in REAL TIME from start of first LFP file; includes gaps between merged files
%fromlfpfilestart is in ms in RELATIVE TIME to start of first LFP file; does not include file gaps. 
% For use in neuroscope

if isfile([basename '.lfptimestamps.mat']) % '_1250Hz_lfptimestamps.mat'
    load([basename '.lfptimestamps.mat'])
    disp(['Loading ' basename '.lfptimestamps.mat'])
else
    filesamplerate = rawlfpsamprate;
    
    if ~isempty(baseLFPfolders)
    clear fID dataArray
    formatSpec = '%q%[^\n\r]';
    % '/home/wahlberg/Exp_Data/M1_Nov2021/20211123/2021-11-23_15-34-56/RecordNode108',...
        i = 1;
        for folder = 1:length(baseLFPfolders)
            d = dir(fullfile(baseLFPfolders{folder}, 'settings*'));

            fID = fopen(fullfile(baseLFPfolders{folder},d(1).name),'r');
            dataArray = textscan(fID,formatSpec,'Delimiter', ',', 'ReturnOnError', false);

            tmp = split(dataArray{1}{6});
            tmp = split(tmp{4},'<');
            lfptimestamps.starttime{i} = tmp{1};
            i = i+1;

            if length(d) > 1 % if two experiments
                clear fID dataArray
                fID = fopen(fullfile(baseLFPfolders{folder},d(2).name),'r');
                dataArray = textscan(fID, formatSpec, 'Delimiter', ',',  'ReturnOnError', false);

                tmp = split(dataArray{1}{6});
                tmp = split(tmp{4},'<');
                lfptimestamps.starttime{i} = tmp{1};
                i = i+1;
            end
            clear fID d dataArray
        end
    end

    if TimeIntervalCombinedpath
        % TimeIntervalCombined is file from utku's merging files code
        formatSpec = '%q%q%q%[^\n\r]';
        fID = fopen(TimeIntervalCombinedpath,'r');
        dataArray = textscan(fID, formatSpec, 'Delimiter', ',',  'ReturnOnError', false);

        filesamplerate = str2double(dataArray{3}{2});
        
        n = 1;
         lfptimestamps.fromlfpclockstart = timetable(seconds(0));

        for i = 2:length(dataArray{2}) % 4 because not using some of the earlier empty files
            lfptimestamps.nsmps(n) = str2double(dataArray{2}{i});
            lfptimestamps.startdatetimes{n} = datetime(dataArray{1}{i},...
                'inputFormat','yyyy-MM-dd HH:mm:ss.SSS','Format','yyyy-MM-dd HH:mm:ss.SSS');
            
            s = linspace(0,lfptimestamps.nsmps(n)/filesamplerate,lfptimestamps.nsmps(n));
            s = seconds(s);

            lfptimestamps.clocktimes{n} = s + lfptimestamps.startdatetimes{n}; % actual time of day
            lfptimestamps.fromclock{n} = lfptimestamps.clocktimes{n} - lfptimestamps.startdatetimes{1}; % time relative to start of first lfpfile

            fromclockvec = timetable((lfptimestamps.fromclock{n})');
            lfptimestamps.fromlfpclockstart = [lfptimestamps.fromlfpclockstart;fromclockvec];
            
            n = n+1;
        end


        %nlfp = length(lfptimestamps.nsmps); % total number of lfp files
        
        % in ms from start of file, not counting for merge gaps
        ntotalsamps=sum(lfptimestamps.nsmps(:));
        lfptimestamps.fromlfpfilestart = (linspace(0,ntotalsamps/filesamplerate*1000,ntotalsamps))';
        lfptimestamps.fromlfpfilestart = [0;lfptimestamps.fromlfpfilestart];
        % in ms (real clock time) from start of first lfp file, accounting for merge gaps
        lfptimestamps.fromlfpclockstart = seconds(lfptimestamps.fromlfpclockstart.Time(:))*1000;
        % NOT getting rid of the original 0 allocation because the
        % mergedLFP file is 1 sample more than the sum of the parts
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

    behaviortimestamps.fromlfpclockstart = timetable(seconds(0));
    behaviortimestamps.fromlfpfilestart = [];
    for f = 1:length(behdir)
        fID = fopen(fullfile(behaviorpath, behdir(f).name),'r');
        dataArray = textscan(fID, formatSpec, 'Delimiter', ',',  'ReturnOnError', false);

        behaviortimestamps.startdatetimes{f} = datetime(dataArray{12}{1},...
            'inputFormat','yyyy-MM-dd hh.mm.ss.SSS a','Format','yyyy-MM-dd HH:mm:ss.SSS');
        s = (seconds(str2double(dataArray{2}(header_length:end))))';
        

        behaviortimestamps.clocktimes{f} = s + behaviortimestamps.startdatetimes{f}; 
        behaviortimestamps.fromclock{f} = behaviortimestamps.clocktimes{f}-lfptimestamps.startdatetimes{1}; % in minutes
        behaviortimestamps.nsmps(f) = str2double(dataArray{16}{1}); % total frames

        fromclockvec = timetable((behaviortimestamps.fromclock{f})');
        behaviortimestamps.fromlfpclockstart = [behaviortimestamps.fromlfpclockstart;fromclockvec];

        behaviortimestamps.framerate{f} = dataArray{8}{1};

        clear dataArray
    end
 filesamplerate = str2double(behaviortimestamps.framerate{1});
    % in ms (real clock time) from start of first lfp file, accounting for merge gaps
       behaviortimestamps.fromlfpclockstart = seconds(behaviortimestamps.fromlfpclockstart.Time(:))*1000;
    behaviortimestamps.fromlfpclockstart(1) = []; % to get rid of the original zero allocation
   
% Align with lfp

% works but way too slow, 35 min
disp('Starting creation of behaviortimestamps.indexintolfp; will take awhile')
behaviortimestamps.indexintolfp = zeros(size(behaviortimestamps.fromlfpclockstart));
for b= 1:length(behaviortimestamps.fromlfpclockstart)
    tic
    difvec = lfptimestamps.fromlfpclockstart - behaviortimestamps.fromlfpclockstart(b);
    [~,index] = min(abs(difvec));
    behaviortimestamps.indexintolfp(b) = index;
    T(b) = toc;
end

disp(['Time taken to create behaviortimestamps.indexintolfp = ' sum(tT)/50 ' minutes'])
   
behaviortimestamps.fromlfpfilestart = lfptimestamps.fromlfpfilestart(behaviortimestamps.indexintolfp);
end

% lfptimestamps.notes.alignment = 'time 0 is the starttime of the first LFP file';
% lfptimestamps.notes.samplingrate = lfpsamplingrate;
% 
% behaviortimestamps.notes.alignment = 'time 0 is the starttime of the first LFP file';
% behaviortimestamps.notes.samplingrate = behaviorsamplingrate;

% include notes about conversion, original sampling rates, etc

%% Make field where 1s are the good time periods and 0 are the bad time periods

remove = removetimestamps.sampformat;

lfptimestamps.keeptimestamps = ones(size(lfptimestamps.fromlfpclockstart));
for r = 1:length(remove)
    lfptimestamps.keeptimestamps(remove(r,1):remove(r,2)) = 0;
end

behaviortimestamps.keeptimestamps = ones(size(behaviortimestamps.fromlfpclockstart));
for r = 1:length(behaviortimestamps.indexintolfp)
    if isnan(lfptimestamps.fromlfpfilestart(behaviortimestamps.indexintolfp(r)))
        behaviortimestamps.keeptimestamps(r) = 0;
    end
end

%% Check that it looks right! Also you can visually check if any files seem to be missing

if showplot == 1

    figure
    title('Check from LFP clock start')
    hold all; ylim([-1 2]); xlabel('ms from LFP clock start')

    xvec = ones(length(find(behaviortimestamps.keeptimestamps == 1)),1)+0.5;
    scatter((behaviortimestamps.fromlfpclockstart(behaviortimestamps.keeptimestamps == 1)'),xvec',0.5,'g')

    xvec = ones(length(behaviortimestamps.fromlfpclockstart),1);
    scatter(behaviortimestamps.fromlfpclockstart',xvec',0.5,'b')

    xvec = zeros(length(lfptimestamps.fromlfpclockstart),1);
    scatter(lfptimestamps.fromlfpclockstart',xvec',0.5,'r')

    xvec = zeros(length(find(lfptimestamps.keeptimestamps == 1)),1)-0.5;
    scatter((lfptimestamps.fromlfpclockstart(lfptimestamps.keeptimestamps == 1)')',(xvec)',0.5,'k')
  
    legend({'Beh ts to use','All beh ts','All LFP ts','LFP ts to use'})
    
    hold off
    figure
    title('Check from LFP file start')
    hold all; ylim([-1 2]); xlabel('ms from LFP file start')

    xvec = ones(length(behaviortimestamps.fromlfpclockstart),1);
    scatter(lfptimestamps.fromlfpfilestart(behaviortimestamps.indexintolfp)',xvec',0.5,'b')
    
    xvec = zeros(length(lfptimestamps.fromlfpfilestart),1);
    scatter(lfptimestamps.fromlfpfilestart',xvec',0.5,'r')

end

%% Save

if savelfpfilename
    save(savelfpfilename,'lfptimestamps')
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




