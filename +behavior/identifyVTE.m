function exporttrials = identifyVTE(varargin)

p = inputParser;
addRequired(p,'trials',@isstruct)
addRequired(p,'boundaries',@isstruct)
addParameter(p,'updatetrials',false,@islogical)
addParameter(p,'checkplot',false,@islogical)

parse(p,varargin{:})
trials = p.Results.trials;
boundaries = p.Results.boundaries;
updatetrials = p.Results.updatetrials;
checkplot = p.Results.checkplot;

tmptrials = trials;

%% Initial vals
nTrials = size(tmptrials,1);
bottomchoicebox = boundaries.choice.x(1);  % same as boundary for goal_exit
choice_exitL = boundaries.goal_l.y(1);
choice_exitR = boundaries.goal_r.y(2);

%% Identify when the animal is in the choice box

for tri = 1:nTrials
    x = tmptrials.xtrialdata{tri};
    y = tmptrials.ytrialdata{tri};

    % when animal enters choice box (determines firstpass of defining trials)
    xind = find(x < bottomchoicebox & y > boundaries.center.y(1) & y < boundaries.center.y(2)); % when animal is in box
    if xind(1) == 1; xind(1) = []; end
    if xind(end) == length(y); xind(end)=[]; end
    xpost = xind + 1; % index for time right before in box
    tmp = find(x(xpost) > bottomchoicebox);
    choiceenter = xpost(tmp)-1;

    if strcmp(tmptrials.trialDirection{tri},'Left') == 1

        % when animal exits choice box (LEFT)
        yind = find(x > bottomchoicebox & ...
            y > boundaries.center.y(1) & y < boundaries.center.y(2)); % animal is in box
        if yind(1) == 1; yind(1) = []; end
        if yind(end) == length(x); yind(end)=[]; end
        ypost = yind + 1;
        tmp = find(y(ypost)>boundaries.center.y(2));
        choiceexit = ypost(tmp)-1;

        %%%%% Identify VTE %%%%%%%%%%

        ychoice = y(choiceenter:choiceexit); % only need y
        ydiff = diff(ychoice);
        if ydiff < 0 ; tmptrials.vte{tri} = 1;
        else; tmptrials.vte{tri} = 0; end

    else % Rightward trial

        % when animal exits choice box (RIGHT)
        yind = find(x > bottomchoicebox & ...
            y > boundaries.center.y(1) & y < boundaries.center.y(2)); % animal is in box
        if yind(1) == 1; yind(1) = []; end
        if yind(end) == length(x); yind(end)=[]; end
        ypost = yind + 1;
        tmp = find(y(ypost)<boundaries.center.y(1));
        choiceexit = ypost(tmp)-1;

        %%%%% Identify VTE %%%%%%%%%%

        ychoice = y(choiceenter:choiceexit); % only need y
        ydiff = diff(ychoice);
        if ydiff > 0 ; tmptrials.vte{tri} = 1;
        else; tmptrials.vte{tri} = 0; end
    end
end

if checkplot == 1
    figure
    for tri = 1:nTrials
        plot(trials.xtrialdata{tri},trials.ytrialdata{tri})
        title(['VTE: ' num2str(tmptrials.vte{tri})])
        pause
    end
end

if updatetrials == 1; exporttrials = tmptrials;
else; exporttrials = trials; end




