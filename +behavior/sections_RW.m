function bounds = sections_RW(varargin)
% [bounds, rot_x, rot_y, rot_ang] = sections(x, y, skip_rot_check, ...)
%
% Will Mau, Nat Kinsky
% Edited Rachel Wahlberg Mar 2022
%   This function takes position data and partitions the maze into
%   sections.
%
%   INPUTS:
%       X and Y: Position vectors
%
%   OUTPUTS:
%       BOUNDS: Struct containing coordinates for the corners of maze
%       sections in the following intuitively-named fields.
%           base = Start position (vertical stripes side).
%           center = Middle stem.
%           choice = Choice point (triangle side).
%           approach_l = Approaching left arm.
%           approach_r = Approaching right arm.
%           left = Left arm.
%           right = Right arm.
%           return_l = Returning to start position from left arm.
%           return_r = Returning to start position right right arm.
%           goal_l = in reward zone on left arm.
%           goal_r = in reward zone on right arm.
%

p = inputParser;
addRequired(p,'x',@isnumeric);
addRequired(p,'y',@isnumeric);
addParameter(p,'boundariesplotcheck',true,@islogical);
addParameter(p,'savefigpath',[],@ischar);

parse(p,varargin{:});
x = p.Results.x;
y = p.Results.y;
boundariesplotcheck = p.Results.boundariesplotcheck;
savefigpath = p.Results.savefigpath;

%% Get xy coordinate bounds for maze sections.
ymax = max(y); xmin = min(x);
xmax = max(x); ymin = min(y);

%%
rest_x = 1/4*range(x) + min(x);% to min(y) for bottom of box
choice_enter = 3/4*range(x) + min(x); % to max(y) for top of box
choice_exitR = 1/3*range(y) + min(y);
choice_exitL = 2/3*range(y) + min(y);
center_R = 2/5*range(y) + min(y);
center_L = 3/5*range(y) + min(y);
goal_R = 1/4*range(y) + min(y);
goal_L = 3/4*range(y) + min(y);

%% Define bounds

% y vals are BOTTOM TOP
% x vals are LEFT RIGHT
center.y = [center_R,center_L,]; % why are they repeated again?
center.x = [rest_x,choice_enter];

% left arm (between choice point and goal)
approach_l.y = [choice_exitL,goal_L];
approach_l.x = [choice_enter,xmax];

% right arm (between choice point and goal)
approach_r.y = [goal_R,choice_exitR];
approach_r.x = [choice_enter,xmax];

%left return
return_l.y = [center_L,ymax];
return_l.x = [rest_x,choice_enter];

%Right return
return_r.y = [ymin,center_R];
return_r.x = [rest_x,choice_enter];

%left goal
goal_l.y = [goal_L,ymax];
goal_l.x = [choice_enter,xmax];

%right goal
goal_r.y = [ymin,goal_R];
goal_r.x = [choice_enter,xmax];

%Choice
choice.y = [choice_exitR,choice_exitL];
choice.x = [choice_enter,xmax];

%Rest box
rest.y = [ymin,ymax];
rest.x = [xmin,rest_x];

%Centroid of maze, for linearization purposes
centroid.y = mean(ymin+ymax);
centroid.x = mean(xmin+xmax);

%% Check with plot

if boundariesplotcheck == 1
    plotthese = {center,choice,approach_l,approach_r,return_l,return_r,goal_l,goal_r,rest};
    colors = {'k','k','g','g','b','b','y','y','r'};

    figure
    hold all
    plot(x,y)
    for i = 1:9
        var = plotthese{i};
        line([var.x(1),var.x(2)],[var.y(1),var.y(1)],'Color',colors{i})
        line([var.x(1),var.x(2)],[var.y(2),var.y(2)],'Color',colors{i})
        line([var.x(1),var.x(1)],[var.y(1),var.y(2)],'Color',colors{i})
        line([var.x(2),var.x(2)],[var.y(1),var.y(2)],'Color',colors{i})
        %pause
    end

    if ~isempty(savefigpath)
        save(gcf,savefigpath)
    end


end

%% Output.
bounds.center = center;
bounds.choice = choice;
bounds.approach_l = approach_l;
bounds.approach_r = approach_r;
bounds.return_l = return_l;
bounds.return_r = return_r;
bounds.goal_l = goal_l;
bounds.goal_r = goal_r;
bounds.rest = rest;
bounds.centroid = centroid;
bounds.range.y = [ymin ymax];
bounds.range.x = [xmin xmax];
end
