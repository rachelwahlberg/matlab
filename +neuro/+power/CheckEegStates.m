function gCheckEegStates = CheckEegStates(varargin) % fileinfo was a previous arg
%function CheckEegStates(fileinfo)

p = inputParser;
addRequired(p,'basepath',@ischar)
addRequired(p,'basename',@ischar)
addParameter(p,'lfp',[],@isnumeric) % if you want to send in data directly instead of loading
%addParameter(p,'lfpCh',[],@isnumeric)
addParameter(p,'redo_flag',0,@islogical)
addParameter(p,'State',[],@istruct) % I think it's a struct?
addParameter(p,'AuxData',[],@isnumeric) % I think it's a matrix? 
addParameter(p,'MinLen',5,@isnumeric) % in seconds
addParameter(p,'Window',1,@isnumeric), % is in seconds - for spectrogram calculatuions
addParameter(p,'FreqRange',[1 100], @isnumeric) % Hz - range of freq for spectrogram

parse(p,varargin{:})
basepath = p.Results.basepath;
basename = p.Results.basename;
lfp = p.Results.lfp;
%lfpCh = p.Results.lfpCh;
redo_flag = p.Results.redo_flag;
State = p.Results.State;
AuxData = p.Results.AuxData;
MinLen = p.Results.MinLen;
Window = p.Results.Window;
FreqRange = p.Results.FreqRange;

%[redo_flag,State,AuxData] = DefaultArgs(varargin,{0,[],[]});

% auxil. struct for gui
global gCheckEegStates
gCheckEegStates = struct;

cd(basepath);
%basepath = pwd; 
%filename = fileinfo.name;

%basename = [basepath '/' filename]; % '/' filename];

%constants
UseRulesBefore = 0; % flag to switch heuristics for periods cleaning after automatic segmentation
%MinLen=5; %seconds
Par = LoadPar([basename '.xml']);
if isfield(Par,'lfpSampleRate')
    eSampleRate = Par.lfpSampleRate;
else
    eSampleRate = 1250;
end
nChannels = Par.nChannels;
%Window = 1; %sec - for spectrogram calculation
%FreqRange = [1 100]; % Hz - range of freq. for the spectrogram

if ~isempty(State)
    % load segmentation results
    if FileExists([basename '.sts.' State])
        Per = load([basename '.sts.' State]);
    elseif FileExists([basename '.states.res'])
        Per = SelectStates(basename,State,eSampleRate*MinLen);
    else
        Per = [];
    end       
    if isempty(Per)
        fprintf('No segments detected in %s state\n',State)
        %return;
        Per = [];
    end
else
    Per = [];
    State = '';
end

Channels = [];

% if ~isempty(fileinfo.CA1thetach)
%     ca1 = find(fileinfo.CA==1); 
% %     Channels = [Channels fileinfo.GammaCh(ca1(1))];
%     Channels = fileinfo.CA1thetach + 1;
% elseif ~isempty(fileinfo.CA3thetach)
%     ca3 = find(fileinfo.CA==3);
% %     Channels = [Channels fileinfo.GammaCh(ca3(1))];
%     Channels = fileinfo.CA3thetach+1;
%     fprintf('Warning: Using CA3 theta Channel!!\n')
% else
%     warning('No theta channel available!!');
% end

if FileExists([basename '.eegseg.mat']) & ~redo_flag
    load([basename '.eegseg.mat']); % load whitened spectrograms from EegSegmentation results
else
    if isempty(Channels)
    ch = inputdlg({'Enter channels to use'},'Channels',1,{'1'});
    %Channels = num2str(ch{1});
    Channels = str2double(ch{1});
    end
 
    % now compute the spectrogram
    if ~isempty(lfp) % from command line, in matrix nsamps x nch 
        Eeg = lfp(:,Channels);
    else % if isempty(lfp), then load from file. 
        if exist([basename '.eeg'],'file')
            Eeg = readmulti([basename '.eeg'], nChannels, Channels);
        elseif exist([basename '.lfp'],'file')
            Eeg = readmulti([basename '.lfp'], nChannels, Channels);
        else
            keyboard
        end
    end

% highband = 600; % bandpass filter range
% lowband = 15; % (250Hz to 80Hz)
% 
%     forder = 100;  % filter order has to be even. The longer the more selective, but the operation
%                % will be linearly slower to the filter order. 100-125 for 1.25Hz (500 for 5 KH
%     % avgfilorder = 101; % should not change this. length of averaging filter
%     forder = ceil(forder/2)*2; % to make sure filter order is even
%     firfiltb = fir1(forder,[lowband/eSampleRate*2,highband/eSampleRate*2]); % calculate convolution func
%     
%     Eeg = Filter0(firfiltb,Eeg); % filtering

    nFFT = 2^round(log2(2^11)); %compute nFFT according to different sampling rates
    SpecWindow = 2^round(log2(Window*eSampleRate));% choose window length as power of two
%     weeg = WhitenSignal(Eeg,[],1,[1.0306 -0.0451]);
    weeg = WhitenSignal(Eeg,100,1,[]);
%    weeg = WhitenSignal(Eeg,[],1,[]);
    %%%%      mtcsglong(x,nFFT,Fs,WinLength,nOverlap,NW,Detrend,nTapers,FreqRange);
    [Pxx,f,t]=neuro.power.mtcsglong(weeg,nFFT,eSampleRate,SpecWindow,[],2,'linear',[],FreqRange);
%     [Pxx,f,t]=mtptcsg(weeg,[],[],nFFT,eSampleRate,SpecWindow,[],2,'linear',[],FreqRange);
    ifsave = 1; %input('Do you want to save the spectrum? [1/0]')
    if ifsave
        save([basename '.eegseg.mat'],'Pxx','f','t','-v6');
    end
end

t = (t(2)-t(1))/2 +t;
   
% computer the/del ratio and detect transitions automatically - not used at
% the momnet, maybe later
%[thratio] = TDRatioAuto(Pxx,f,t,MinLen);
%[thratio, ThePeriods] = TDRatioAuto(Pxx,f,t,MinLen);

%now apply the rules to filter out junk states or make continuous periods
% to be implemented later
if UseRulesBefore
    switch State
        case 'REM'

    end
end

nAuxData = size(AuxData,1);

% fill the global structure for future use
gCheckEegStates.Channels = Channels;
gCheckEegStates.nChannels = nChannels;
gCheckEegStates.FileBase = [basename '_1250Hz'];
gCheckEegStates.State = State;
gCheckEegStates.t = 10; %is seconds
gCheckEegStates.eFs = eSampleRate;
gCheckEegStates.trange = [t(1) t(end)];
gCheckEegStates.Periods = Per/eSampleRate; % in seconds
gCheckEegStates.Mode = 't';
gCheckEegStates.nPlots=length(Channels)+1+nAuxData;%number of plots is fixed to 3: 2 spectrograms and one zoomin trace. If changed reuiqres fixes in _aux
gCheckEegStates.lh=cell(gCheckEegStates.nPlots,1);
gCheckEegStates.Window = Window*eSampleRate*2;
gCheckEegStates.SelLine=[];
gCheckEegStates.cposh=cell(gCheckEegStates.nPlots,1);
gCheckEegStates.FreqRange = FreqRange;
gCheckEegStates.newl=[];
gCheckEegStates.tstep = t(2)-t(1);
gCheckEegStates.coolln = [];
gCheckEegStates.LastBut = 'normal';
gCheckEegStates.nAuxData = nAuxData;
% create and configure the figure
gCheckEegStates.figh = figure('ToolBar','none');
set(gCheckEegStates.figh, 'Position', [3 128 1276 620]); %change Postion of figure if you have low resolution
set(gCheckEegStates.figh, 'Name', 'CheckEegStates : traces');
set(gCheckEegStates.figh, 'NumberTitle', 'off');


% put the uitoolbar and uimenu definitions here .. may require rewriting
% some callbacks as actions rather then cases of actions (e.g. key
% pressing)

%% now do the plots

for ii=1:length(gCheckEegStates.Channels)
    subplot(gCheckEegStates.nPlots,1,ii);
    imagesc(t,f,log(squeeze(Pxx(:,:,ii)))');axis xy; ylim([1 20]);
    hold on
    if ii==1
        title('Spectrogram'); ylabel('Frequency (Hz)');
    end
end

if nAuxData>0
    for ii=[1:nAuxData]
        subplot(gCheckEegStates.nPlots,1,ii+length(gCheckEegStates.Channels));
        DisplayAuxData(AuxData(:,ii));
        hold on
        
    end
end

%  subplot(gCheckEegStates.nPlots,1,2)
%  plot(t,thratio);axis tight;
%  set(gca,'YTick',[]);
%  hold on
%  ylabel('Theta/Delta raio'); xlabel('Seconds');

subplot(gCheckEegStates.nPlots,1,gCheckEegStates.nPlots)
CheckEegStates_aux('traces'); % plot the current eeg traces
hold on

%plots lines
if ~isempty(Per)
CheckEegStates_aux('lines');
end
% assign functions for mouse and keyboard click
set(gCheckEegStates.figh,'WindowButtonDownFcn','CheckEegStates_aux(''mouseclick'')');
set(gCheckEegStates.figh,'KeyPressFcn', 'CheckEegStates_aux(''keyboard'')');

%msave([FileBase '.states.' State],round(ThePeriods*eSampleRate));
return

%% Define functions

function [thratio, ThePeriods] = TDRatioAuto(Pxx,f,t,MinLen)

%automatic theta periods detection just using the thetaratio
thfin = find(f>6 & f<9);
thfout = find(f<5 | (f> 12& f<25));
thratio = log(mean(squeeze(Pxx(:,thfin,1)),2))-log(mean(squeeze(Pxx(:,thfout,1)),2));

if nargout>1
    nStates =2;
    % fit gaussian mixture and to HMM - experimental version .. uses only thetaratio
    [TheState thhmm thdec] = gausshmm(thratio,nStates,1,0);

    for i=1:nStates
        thratio_st(i) = mean(thratio(TheState==i));
    end

    [dummy TheInd] = max(thratio_st);
    InTh = (TheState==TheInd);
    DeltaT = t(2)-t(1);
    MinSeg = round(MinLen/DeltaT);

    TransTime = ThreshCross(InTh,0.5,MinSeg);
    ThePeriods = t(TransTime);
end
return

function DisplayAuxData(Data)

        nEl =  size(Data,2);
        if nEl<4 
            err=1;
        elseif ~isstr(Data{4})
            err=1;
        else
            err=0;
        end
        if err 
            warning('AuxData has to be cell array where each row is : xaxis, yaxis, data, display_func');
            close 
            return;
        end
        
        switch Data{4} %switch by functions
            case 'plot'
                plot(Data{1}, Data{3});
            case 'imagesc'
                if length(Data{1})~=size(Data{3},1) & length(Data{1})~=size(Data{3},2)
                    Data{3}=Data{3}';
                end
                    
                imagesc(Data{1},Data{2}, Data{3}');
            otherwise 
                error('wrong data display function');
        end
        


return