function bounds = sections_RW(x, y, plotcheck)
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

%% Get xy coordinate bounds for maze sections.
xmax = max(x); xmin = min(x);
ymax = max(y); ymin = min(y);

%%
rest_y = 1/4*range(y) + min(y);% to min(y) for bottom of box
choice_enter = 3/4*range(y) + min(y); % to max(y) for top of box
choice_exitL = 1/3*range(x) + min(x);
choice_exitR = 2/3*range(x) + min(x);
center_L = 2/5*range(x) + min(x);
center_R = 3/5*range(x) + min(x);
goal_L = 1/4*range(x) + min(x);
goal_R = 3/4*range(x) + min(x);





%% Define bounds

% y vals are BOTTOM TOP
% x vals are LEFT RIGHT
center.x = [center_L,center_R]; % why are they repeated again?
center.y = [rest_y,choice_enter];

% left arm (between choice point and goal)
approach_l.x = [goal_L,choice_exitL];
approach_l.y = [choice_enter,ymax];

% right arm (between choice point and goal)
approach_r.x = [choice_exitR,goal_R];
approach_r.y = [choice_enter,ymax];

%left return
return_l.x = [xmin,center_L];
return_l.y = [rest_y,choice_enter];

%Right return
return_r.x = [center_R,xmax];
return_r.y = [rest_y,choice_enter];

%left goal
goal_l.x = [xmin,goal_L];
goal_l.y = [choice_enter,ymax];

%right goal
goal_r.x = [goal_R,xmax];
goal_r.y = [choice_enter,ymax];

%Choice
choice.x = [choice_exitL,choice_exitR];
choice.y = [choice_enter,ymax];

%Rest box
rest.x = [xmin,xmax];
rest.y = [ymin,rest_y];

%Centroid of maze, for linearization purposes
centroid.x = mean(xmin+xmax);
centroid.y = mean(ymin+ymax);

%% Check with plot

if plotcheck == 1
    figure
    hold all
    plot(x,y)
    %center box
    line([center.x(1),center.x(1)],[center.y(1),center.y(2)],'Color','k')
    line([center.x(2),center.x(2)],[center.y(1),center.y(2)],'Color','k')
    line([center.x(1),center.x(2)],[center.y(1),center.y(1)],'Color','k')
    line([center.x(1),center.x(2)],[center.y(2),center.y(2)],'Color','k')
    %goal box
    line([choice.x(1),choice.x(1)],[choice.y(1),choice.y(2)],'Color','g')
    line([choice.x(2),choice.x(2)],[choice.y(1),choice.y(2)],'Color','g')
    line([choice.x(1),choice.x(2)],[choice.y(1),choice.y(1)],'Color','g')
    line([choice.x(1),choice.x(2)],[choice.y(2),choice.y(2)],'Color','g')
    %approach left
    line([approach_l.x(1),approach_l.x(1)],[approach_l.y(1),approach_l.y(2)],'Color','b')
    line([approach_l.x(2),approach_l.x(2)],[approach_l.y(1),approach_l.y(2)],'Color','b')
    line([approach_l.x(1),approach_l.x(2)],[approach_l.y(1),approach_l.y(1)],'Color','b')
    line([approach_l.x(1),approach_l.x(2)],[approach_l.y(2),approach_l.y(2)],'Color','b')
    %approach right
    line([approach_r.x(1),approach_r.x(1)],[approach_r.y(1),approach_r.y(2)],'Color','b')
    line([approach_r.x(2),approach_r.x(2)],[approach_r.y(1),approach_r.y(2)],'Color','b')
    line([approach_r.x(1),approach_r.x(2)],[approach_r.y(1),approach_r.y(1)],'Color','b')
    line([approach_r.x(1),approach_r.x(2)],[approach_r.y(2),approach_r.y(2)],'Color','b')
    %return left
    line([return_l.x(1),return_l.x(1)],[return_l.y(1),return_l.y(2)],'Color','y')
    line([return_l.x(2),return_l.x(2)],[return_l.y(1),return_l.y(2)],'Color','y')
    line([return_l.x(1),return_l.x(2)],[return_l.y(1),return_l.y(1)],'Color','y')
    line([return_l.x(1),return_l.x(2)],[return_l.y(2),return_l.y(2)],'Color','y')
    %return right
    line([return_r.x(1),return_r.x(1)],[return_r.y(1),return_r.y(2)],'Color','y')
    line([return_r.x(2),return_r.x(2)],[return_r.y(1),return_r.y(2)],'Color','y')
    line([return_r.x(1),return_r.x(2)],[return_r.y(1),return_r.y(1)],'Color','y')
    line([return_r.x(1),return_r.x(2)],[return_r.y(2),return_r.y(2)],'Color','y')
    %goal left
    line([goal_l.x(1),goal_l.x(1)],[goal_l.y(1),goal_l.y(2)],'Color','r')
    line([goal_l.x(2),goal_l.x(2)],[goal_l.y(1),goal_l.y(2)],'Color','r')
    line([goal_l.x(1),goal_l.x(2)],[goal_l.y(1),goal_l.y(1)],'Color','r')
    line([goal_l.x(1),goal_l.x(2)],[goal_l.y(2),goal_l.y(2)],'Color','r')
    %goal right
    line([goal_r.x(1),goal_r.x(1)],[goal_r.y(1),goal_r.y(2)],'Color','r')
    line([goal_r.x(2),goal_r.x(2)],[goal_r.y(1),goal_r.y(2)],'Color','r')
    line([goal_r.x(1),goal_r.x(2)],[goal_r.y(1),goal_r.y(1)],'Color','r')
    line([goal_r.x(1),goal_r.x(2)],[goal_r.y(2),goal_r.y(2)],'Color','r')
    %rest box
    line([rest.x(1),rest.x(1)],[rest.y(1),rest.y(2)],'Color','k')
    line([rest.x(2),rest.x(2)],[rest.y(1),rest.y(2)],'Color','k')
    line([rest.x(1),rest.x(2)],[rest.y(1),rest.y(1)],'Color','k')
    line([rest.x(1),rest.x(2)],[rest.y(2),rest.y(2)],'Color','k')
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
end
