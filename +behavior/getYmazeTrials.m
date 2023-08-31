 function trials = getYmazeTrials(varargin)
% Maze should be facing vertical with rest area aligned right above y = 0

p = inputParser;
addRequired(p,'boundaries',@isstruct)
addRequired(p,'maze',@isstruct)
addRequired(p,'behaviorsamprate',@isnumeric) % needed if saving an event file!
addParameter(p,'behaviortimestamps',[],@isstruct) % needed if removing bad periods!
addParameter(p,'checkplot',false,@islogical)
addParameter(p,'newtrialfilename',[],@ischar)

parse(p,varargin{:})
boundaries = p.Results.boundaries;
maze = p.Results.maze;
behaviorsamprate = p.Results.behaviorsamprate;
behaviortimestamps = p.Results.behaviortimestamps;
checkplot = p.Results.checkplot;
newtrialfilename = p.Results.newtrialfilename;

%% initial vals
% using interpolated values otherwise the nans mess things up

% if maze isn't rotated,this could lead to errors in enter/exit definitions
if range(maze.position.interpolatedx) > range(maze.position.interpolatedy)
    warning('Did you rotate/shift your maze?')
end

trials = [];

%% Get transitions between boundaries
x = maze.position.interpolatedx;
y = maze.position.interpolatedy;

% get relevant boundaries
bottomchoicebox = boundaries.choice.x(1);  % same as boundary for goal_exit
toprestbox= boundaries.rest.x(2);
goal_enterL = boundaries.goal_l.y(1);
goal_enterR = boundaries.goal_r.y(2);

% when animal enters rest box (start of PRE TRIAL, end of POST TRIAL)
xind = find(x > toprestbox); % & x > boundaries.center.x(1) & x < boundaries.center.x(2)); % when animal has left rest box
if xind(1) == 1; xind(1) = []; end
if xind(end) == length(y); xind(end)=[]; end
xpost = xind + 1; % index for time right before in box
tmp = find(x(xpost) < toprestbox);
restenter = xpost(tmp)-1;

% when animal exits rest box (end of PRETRIAL, start of TRIAL)
xind = find(x < toprestbox); % & x > boundaries.center.x(1) & x < boundaries.center.x(2)); % when animal has left rest box
if xind(1) == 1; xind(1) = []; end
if xind(end) == length(y); xind(end)=[]; end
xpost = xind + 1; % index for time right before in box
tmp = find((xpost) > toprestbox);
restexit = xpost(tmp)-1;

% when animal enters choice box (determines firstpass of defining trials)
xind = find(x < bottomchoicebox & y > boundaries.center.y(1) & y < boundaries.center.y(2)); % when animal is in box
if xind(1) == 1; xind(1) = []; end
if xind(end) == length(y); xind(end)=[]; end
xpost = xind + 1; % index for time right before in box
tmp = find(x(xpost) > bottomchoicebox);
choiceenter = xpost(tmp)-1;

% exiting left goal box
xindL = find(x > bottomchoicebox & y > goal_enterL); %before exiting goal box
if xindL(end) == length(x); xindL(end) = []; end
if xindL(end) == length(y); xindL(end)=[]; end
xpostL = xindL + 1;
tmp = find(x(xpostL) < bottomchoicebox); %find transition
xgoalexitL = xpostL(tmp)-1;

% exiting right goal box
xindR = find(x > bottomchoicebox & y < goal_enterR); %before exiting goal box
if xindR(end) == length(x); xindR(end) = []; end
if xindR(end) == length(y); xindR(end)=[]; end
xpostR = xindR + 1;
tmp = find(x(xpostR) < bottomchoicebox); %find transition
xgoalexitR = xpostR(tmp)-1;

% all exits of left or right goal boxes (end of TRIAL, start of POST TRIAL)
allgoalexits = unique([xgoalexitL xgoalexitR]);

%% Get trial boundaries

%%%%%%%%%%%%%%%%%%%%%%%%%%%% firstpass %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%firstpass(c,1) = closest time relative to (c,2) that rat leaves rest box
%firstpass(c,2) = unique timestamps for when rat crossed into choice box
%firstpass(c,3) = closest time relative to (c,2) that rat leaves a goal box
nanidx = 1;
for c = 1:length(choiceenter)
    comparisons = choiceenter(c) - allgoalexits;
    comparisons(comparisons > 0) = nan;
    closestgoalexit = find(abs(comparisons) == min(abs(comparisons)));
    firstpass(c,2) =  choiceenter(c);

    if isempty(closestgoalexit)
        firstpass(c,3) = nan;
        cnan(nanidx) = c; nanidx =nanidx+1;
    else
    firstpass(c,3) = allgoalexits(closestgoalexit);
    end

    comparisons = choiceenter(c) - restexit;
    comparisons(comparisons < 0) = nan;
    closestrestexit = find(abs(comparisons) == min(abs(comparisons)));
    if isempty(closestrestexit)
        firstpass(c,1) = nan;
        cnan(nanidx) = c; nanidx =nanidx+1;
    else
    firstpass(c,1) = restexit(closestrestexit);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%% secondpass %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% secondpass deleting rules:
% rest exits same, but goal exits diff - ran twice over choice + goal , delete second set of choice/goal
% rest exits same, goal exits same - ran twice over choice, delete second choice rerun
% rest exits diff, goal exits same - ran over rest exit twice, delete first rest exit
d = 1;
for p = 1:length(firstpass)
    if p == length(firstpass)
        break
    end
    if firstpass(p,1) ~= firstpass(p+1,1) %if two subsequent rest exits are different
        if firstpass(p,3) ~= firstpass(p+1,3) %if the goal exits are different
            continue
        else % if two goal exits are the same
            deletepass(d)=p;
            d = d+1;
        end
    else % if two rest exits are the same
        if firstpass(p,3) ~= firstpass(p+1,3) %if the goal exits are different
            deletepass(d) = p+1;
            d = d+1;
        else % if the goal exits are the same
            deletepass(d) = p + 1;
            d = d + 1;
        end
    end
end

deletethese = [unique(deletepass) cnan];
secondpass = firstpass;
secondpass(deletethese,:) = [];

% Quick check to make sure it worked
sprime = secondpass';
catsprime = [sprime(:)];
secondpassdiff = diff(catsprime);
disp(['All timestamps in order: ' num2str(length(secondpassdiff(:)>0) == length(secondpassdiff))])% make sure that timestamps all get bigger

%% Determine L vs R trials
% Determine whether trials were left or right
for trial = 1:length(secondpass)
    if y(secondpass(trial,3)) < goal_enterR
        trialDirection{trial,1} = 'Right';
    else
        trialDirection{trial,1} = 'Left';
    end
end

%% Get trial position data

for trial = 1:length(secondpass)
    xtrialdata{trial,1} = x(secondpass(trial,1):secondpass(trial,3));
    ytrialdata{trial,1} = y(secondpass(trial,1):secondpass(trial,3));
end

%% Get pretrial and posttrial data
% though they kinda overlap anyway in significance?

%pretrial = when rat enters rest box to exits rest box
%posttrial = when rat leaves goal box to enters rest box
ntrials = size(secondpass,1);
startTrial_index = secondpass(:,1); % REST EXIT
choiceEnter_index= secondpass(:,2);
endTrial_index = secondpass(:,3); % GOAL EXIT

for trial = 1:length(endTrial_index)

    if trial == 1
        startPre_index(trial,1) = 1;
        endPre_index(trial,1) = startTrial_index(1)-1;
    else
        comparisons = endTrial_index(trial-1) - restenter;
        comparisons(comparisons > 0) = nan;
        closestrestenter = find(abs(comparisons) == min(abs(comparisons)));
        startPre_index(trial,1) = restenter(closestrestenter); % END of posttrial, START of pretrial
        endPre_index(trial,1) = startTrial_index(trial)-1;
    end
        xpretrial{trial,1} = x(startPre_index(trial):endPre_index(trial));
        ypretrial{trial,1} = y(startPre_index(trial):endPre_index(trial));
end

for trial = 1:length(endTrial_index)
    if trial == length(secondpass)
        startPost_index(trial,1) = endTrial_index(trial)+1;
        endPost_index(trial,1) = length(x);
    else
        startPost_index(trial,1) = endTrial_index(trial)+1;
        endPost_index(trial,1) = startPre_index(trial+1)-1;
    end
    xposttrial{trial,1} = x(startPost_index(trial):endPost_index(trial));
    yposttrial{trial,1} = y(startPost_index(trial):endPost_index(trial));
end

%%%%%%%%%%%%%%%%%%%%% Determine bad trials %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% CHECK

keep = behaviortimestamps.keeptimestamps;
keep(keep == 0) = nan;

for tri = 1:length(startPre_index)

    preepoch = keep(startPre_index(tri):endPre_index(tri));
    if any(isnan(preepoch))
        CleanPre(tri,1) = 0;
    else
        CleanPre(tri,1) = 1;
    end

    triepoch = keep(startTrial_index(tri):endTrial_index(tri));
    if any(isnan(triepoch))
        CleanTrial(tri,1) = 0;
    else
        CleanTrial(tri,1) = 1;
    end

    postepoch = keep(startPost_index(tri):endPost_index(tri));
    if any(isnan(postepoch))
        CleanPost(tri,1) = 0;
    else
        CleanPost(tri,1) = 1;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%% Create table %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% put into ms clock time relative to LFP
startPre_clock = behaviortimestamps.fromlfpclockstart(startPre_index);
endPre_clock = behaviortimestamps.fromlfpclockstart(endPre_index);
startTrial_clock = behaviortimestamps.fromlfpclockstart(startTrial_index);
choiceEnter_clock = behaviortimestamps.fromlfpclockstart(choiceEnter_index);
endTrial_clock = behaviortimestamps.fromlfpclockstart(endTrial_index);
startPost_clock = behaviortimestamps.fromlfpclockstart(startPost_index);
endPost_clock = behaviortimestamps.fromlfpclockstart(endPost_index);

% also put into time relative to lfpfilestart (for neuroscope, spike rasters, etc)
startPre_file = behaviortimestamps.fromlfpfilestart(startPre_index);
endPre_file = behaviortimestamps.fromlfpfilestart(endPre_index);
startTrial_file = behaviortimestamps.fromlfpfilestart(startTrial_index);
choiceEnter_file = behaviortimestamps.fromlfpfilestart(choiceEnter_index);
endTrial_file = behaviortimestamps.fromlfpfilestart(endTrial_index);
startPost_file = behaviortimestamps.fromlfpfilestart(startPost_index);
endPost_file = behaviortimestamps.fromlfpfilestart(endPost_index);

trials = table(startPre_clock,endPre_clock,startTrial_clock,choiceEnter_clock,...
    endTrial_clock,startPost_clock,endPost_clock,startPre_file,endPre_file,...
    startTrial_file,choiceEnter_file,endTrial_file,startPost_file,endPost_file, ...
    trialDirection,CleanPre,CleanTrial,CleanPost,xpretrial,ypretrial,...
    xtrialdata,ytrialdata,xposttrial,yposttrial);

%% Check figures
if checkplot == 1
    figure
    hold all
    for trial = 1:length(secondpass)
        trial
        plot(trials.trialdata{trial},trials.ytrialdata{trial})
        pause
    end
end

    figure
    hold all
    %for trial = 1:length(secondpass)
        for trial = 1:18
        plot(-trials.ytrialdata{trial},trials.xtrialdata{trial})
    end



%% Save mat 
if newtrialfilename
    if strcmp(newtrialfilename(end-3:end),'.mat') == 1 % saving in seconds relative to LFP! 
        save(newtrialfilename,'trials')
        disp([newtrialfilename ' saved.'])
    elseif strcmp(newtrialfilename(end-3:end),'.evt') == 1 % event file
% saving just trial start, choice point, trial end in format for neuroscope

        finalpass = secondpass'/behaviorsamprate * 1000; % into ms
        events.time = single(vertcat(finalpass(:)));

        events.description = cell(length(events.time),1);
        for i = 1:3:length(events.description)
            events.description{i} = "Start";
            events.description{i+1} = "StartChoice"; %trying alphabetical?
            events.description{i+2} = "Stop";
        end

        SaveEvents_RW(newtrialfilename,events) %a buzcode function
    end
end


