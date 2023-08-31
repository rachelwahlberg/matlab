function [trials,linearizedTrials] = gettrials(varargin)

p = inputParser;
addRequired(p,'maze',@isstruct) % comes in in optitrack format
addRequired(p,'behaviortimestamps',@isstruct)
addRequired(p,'mazetype',@ischar)
addParameter(p,'savetrials',false,@islogical)
addParameter(p,'savetrialtype','evt',@ischar) % event file or mat file
addParameter(p,'manualpath',[],@ischar)
addParameter(p,'getlinearizedtrials',false,@islogical)
addParameter(p,'boundariesplotcheck',true,@islogical)
addParameter(p,'plottrials',false,@islogical)
addParameter(p,'basepath',[],@ischar)
addParameter(p,'basename',[],@ischar)

parse(p,varargin{:})
maze = p.Results.maze; 
behaviortimestamps = p.Results.behaviortimestamps;
mazetype = p.Results.mazetype;
savetrials = p.Results.savetrials;
savetrialtype = p.Results.savetrialtype;
manualpath = p.Results.manualpath;
getlinearizedtrials = p.Results.getlinearizedtrials;
boundariesplotcheck = p.Results.boundariesplotcheck;
plottrials = p.Results.plottrials;
basepath = p.Results.basepath;
basename = p.Results.basename;

%% Alternation Y maze
if mazetype == 'alternationYmaze'
    if manualpath  % would not currently work
        load(manualpath)
    end

    boundaries = sections_RW(maze.position.interpolatedx,maze.position.interpolatedy);
    behaviorsamprate = str2double(behaviortimestamps.framerate{1});

    trials = getYmazeTrials(boundaries,maze,behaviorsamprate, ...
        'behaviortimestamps',behaviortimestamps); %'newtrialfilename',newtrialfilename,

    %% LINEARIZE LATER
    % Linearization isn't working well for this rat M1; pausing too often.
    % it works technically (not straight calling) but doesn't account for jumps;
    % when you return to this separate each separate arm of the maze and
    % linearize for each, and then concatenate
    if getlinearizedtrials == 1
        linearizedTrials = LinearizeYmaze(maze,trials,boundaries,mazetype);
    end

    if savetrials
        if strcmp(savetrialtype,'evt') == 1
            disp('saving trial event file')
            newtrialfilename = fullfile(basepath,[basename '.tri.evt']);
            events = trial2nscformat(trials);
            SaveEvents_RW(newtrialfilename,events); % may need to check ms conversion
        elseif strcmp(savetrialtype,'mat') == 1
            disp('saving trial mat file')
            newtrialfilename = fullfile(basepath,[basename '.trials.mat']);
            save(newtrialfilename,'trials');
        else
            disp('specify evt or mat filetype')
        end
    end

end

if plottrials
    % plot trials
    figure
    hold all
    for tri = 1:length(trials.xtrialdata)
        plot(trials.xtrialdata{tri},trials.ytrialdata{tri})
      %  pause
    end
end

%%% Get trial lengths
% e = trials.endTrial_file;
% s = trials.startTrial_file;
% 
% len=(e-s)/1000/60;
% 
% e = trials.endTrial_clock;
% s = trials.startTrial_clock;
% 
% lenC=(e-s)/1000/60;


