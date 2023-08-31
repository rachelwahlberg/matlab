% For old unused code that probably will never be helpful but who knows

if length(lfp.data) < 10000 % if a smaller file
[cleanlfp,removetimestamps] = glitchdetector(lfp,session.extracellular.srLfp); %,'eventfilename',[basepath basename, '.bad.evt']);
[cleanlfp,badCh] = removebadchannels(cleanlfp);
else
baseline = getbaseline(lfp,session.extracellular.srLFP);
winsize = 30*60*session.extracellular.srLFP; % 30 min window chunks
wstart = 1;

cleanlfp = zeros(size(lfp.data));
removetimestamps = [];
for i = 1:length(lfp.data)
win = wstart:wstart+winsize;
lfpwin = mean(lfp.data(win),2);
[cleanwindow,rmtimesWindow] = glitchdetector_chunks(lfpwin,session.extracellular.srLFP);

cleanlfp(win,:) = cleanwindow;
removetimestamps = [removetimestamps; rmtimesWindow];

wstart = wstart + winsize;
end


%% 

% old code
% not useful but you don't wanna throw away

%% Set the start of the first lfp file as time zero

timezero = lfpfile.starttime{1};
%timezero = str2double(split(lfpfile.starttime{1},':')); 
timeend = str2double(split(lfpfile.starttime{1},':'));

% Get lfp starttimes relative to FIRST LFP start time in ms
for file = 1:length(lfpfile.starttime) 
    start = str2double(split(lfpfile.starttime{file},':'));
    lfp_offset(file) = (start(1)-timezero(1))*3600000+(start(2)-timezero(2))*60000+(start(3)-timezero(3))*1000; %start of lfp file relative to start of FIRST LFP file
end

% Get beh starttimes relative to FIRST LFP start time in ms
for file = 1:length(behavior.starttime)
    start = str2double(split(behavior.starttime{file},'.'));
    beh_offset(file) = (start(1)-timezero(1))*3600000+(start(2)-timezero(2))*60000+(start(3)-timezero(3))*1000+start(4); %adding start(4) correct? 
end

lfp_sr = 1250;
opti_sr = 60;
full_length = lfp_offset(end)+nsmps(end)/lfp_sr*1000; % full length in real time (ms) of beginning of first lfp recording to end of last lfp recording
full_vec = linspace(1,full_length,full_length/1000*lfp_sr);

newtimestamps = nan(1,full_length/1000*lfp_sr); % # of samples if there had been constant data collection at 1250 Hz for the entire realtime
