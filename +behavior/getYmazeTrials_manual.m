function trials = getYmazeTrials(manualpath,video_start,boundaries,behavior_timestamps,lfp_timestamps,maze,showfig)
% Maze should be facing vertical with rest area aligned right above y = 0

% get trial timestamps relative to first lfp file
if usemanualscores == 1 % DON'T EXPECT TO WORK CURRENTLY
manual = readtable(manualpath); %the xlsx version stopped being able to load?
lfp_start = lfp_timestamps.starttime{1};

tstart_raw = manual.Start_definedByEntranceIntoRestPeriodEndOfLastTrial_;
tend_raw = manual.Stop_ifPanelDropped_; %duration format, in format of time passed since video start
tstart = tstart_raw + video_start - lfp_start; %trial starts relative
tend = tend_raw + video_start - lfp_start; %trial starts relative

del = isnan(tstart) | isnan(tend);
tstart(del) = [];
tend(del) = [];

if any(tstart > tend)
    error('Check trial start/end times')
end

ntrials = length(tstart);

xpos = maze.position.x;
ypos = maze.position.y;

x = [];
y = [];
for t = 1:ntrials
    ind = behavior_timestamps.reltime_merged.Time > tstart(t) & behavior_timestamps.reltime_merged.Time < tend(t);
    x{t} = xpos(ind);
    y{t} = ypos(ind);
    clear ind
end

end

%% initial vals

x = maze.position.x;
y = maze.position.y;

% if maze isn't rotated,this could lead to errors in enter/exit definitions 
if range(x) > range(y)
    warning('Did you rotate/shift your maze?')
end

trials = [];


%% Get decision entry points
% within a trial period, get when the animal first enters the box to when they
% first exit. put these into cells. Then run the loop below for each trial, get
% out a phi score for each.

% boundaries of choice box
choice_enter = boundaries.choice.y(1); 
choice_exitL = boundaries.choice.x(1);
choice_exitR = boundaries.choice.x(2);






n = 1; th = 1;
%xdecision = boundaries.xenter;
%ydecision = boundaries.yexit;
for t = 1:ntrials
  
    % when animal enters box
    yind = find(y > xdecision); % when animal is in box
    if yind(1) == 1; yind(1) = []; end
    xpre = yind - 1; % index for time right before in box
    tmp = find(x{t}(xpre) < xdecision);
    xtransition = xpre(tmp)+1;

    % left exits
    xind_l = find(y{t} > ydecision(1) & x{t} > xdecision); %inside box
    if xind_l(end) == length(y{t}); xind_l(end) = []; end
    ypost_l = xind_l + 1;
    tmp = find(y{t}(ypost_l) < ydecision(1)); %find transition
    xtrans_l = ypost_l(tmp)-1;

    % right exits
    xind_r = find(y{t} < ydecision(2) & x{t} > xdecision); % inside box
    if xind_r(end) == length(y{t}); xind_r(end) = []; end
    ypost_r = xind_r + 1;
    tmp = find(y{t}(ypost_r) > ydecision(2));
    xtrans_r = ypost_r(tmp)-1;

    % all y exits
    ytransition = unique([xtrans_l xtrans_r]);

    if isempty(xtransition) | isempty(ytransition)
        trials.thrown(th) = t;
        th = th+1;
        continue
    end

  % stem, in index version for linearization code
    center = find(x{t} > cstart & y{t} ); 
    cexit = find(x{t} < xdecision)


    trials.timestamp(n,1) = xtransition(1); % ENTER
    trials.timestamp(n,2) = ytransition(find(ytransition == min(ytransition))); % EXIT

    if t == 15 % special case, overwrite
        trials.timestamp(n,1) = xtransition(2);
    end

    if ismember(trials.timestamp(n,2),xtrans_l)
        trials.direction{n} = 'Left';
    else
        trials.direction{n} = 'Right';
    end

    trials.xdata{n} = x{t}(trials.timestamp(n,1):trials.timestamp(n,2));
    trials.ydata{n} = y{t}(trials.timestamp(n,1):trials.timestamp(n,2));
    n = n+1;
end

ntrials = length(trials.xdata);
x(trials.thrown) = [];
y(trials.thrown) = [];

% plot it out to confirm it looks right
if showfig == 1
    figure
    hold all
    for t = 1:ntrials
        if isnan(trials.timestamp(t,1))
            continue
        end
        plot(x{t}(trials.timestamp(t,1):trials.timestamp(t,2)),y{t}(trials.timestamp(t,1):trials.timestamp(t,2)))
    end
end
