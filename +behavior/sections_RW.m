function bounds = directionsCircularMaze(positionData)

%Send in PositionData object. 


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

% p = inputParser;
% addRequired(p,'x',@isnumeric);
% addRequired(p,'y',@isnumeric);
% addParameter(p,'boundariesplotcheck',true,@islogical);
% addParameter(p,'savefigpath',[],@ischar);

% parse(p,varargin{:});
% x = p.Results.x;
% y = p.Results.y;
% boundariesplotcheck = p.Results.boundariesplotcheck;
% savefigpath = p.Results.savefigpath;

x = positionData.data.X;
z = positionData.data.Z;

%% Get xy coordinate bounds for maze sections.
zmax = max(z); xmin = min(x);
xmax = max(x); zmin = min(z);
zmid = 0.5*(zmin + zmax);
xmid = 0.5*(xmin + xmax);


%% Check with plot

if boundariesplotcheck == 1
    plotthese = {xmin,xmid,xmax,zmin,zmid,zmax};
    colors = {'k','g','b','y','r','o'};

    figure
    hold all
    plot(x,z)
    for i = 1:6
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
bounds.range.y = [zmin zmax];
bounds.range.x = [xmin xmax];
end
