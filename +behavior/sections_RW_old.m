function [bounds, rot_x, rot_y, rotang] = sections_RW(x, y, skip_rot_check, varargin)
% [bounds, rot_x, rot_y, rot_ang] = sections(x, y, skip_rot_check, ...)
%   
% Will Mau, Nat Kinsky
% Edited Rachel Wahlberg Mar 2022
%   This function takes position data and partitions the maze into
%   sections. 
%
%   INPUTS: 
%       X and Y: Position vectors after passing through
%       PreProcessMousePosition. 
%       skip_rot_check = 0(default if blank): perform manual check of rotation of
%           maze, 1: skip manual check - use this if you know you have already
%           performed a rotation and have a previously saved 'rotated.mat' file
%           in the working directory
%       manual_rot_overwrite(optional): default = 0, 1 will perform manual
%       rotation even if 'rotated.mat' exists in the working directory.
%       Use this to overwrite previously performed rotations. Sample use:
%       sections(x,y,0,'manual_rot_overwrite',1).,
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

%% Get varargins
manual_rot_overwrite = 0; % default value
for j = 1:length(varargin)
   if strcmpi(varargin{j},'manual_rot_overwrite') 
      manual_rot_overwrite = varargin{j+1}; 
   end
end
%% Assign skip_rot_check if not specified
if ~exist('skip_rot_check','var') || ~exist(fullfile(pwd,'rotated.mat'),'file')
    skip_rot_check = 0;
end
%% Correct for rotated maze. 
% skewed = 1;
% while skewed
%     
%     if exist(fullfile(pwd,'Pos_align.mat'),'file') % Skip rotating if already done.
%         [rot_x,rot_y,rotang] = rotate_traj(x,y,0);
%         disp('ASSUMING DATA ALREADY ROTATED: Pos_align.mat file found')
%     else
%         %Try loading previous rotation angle.
%         try
%             load(fullfile(pwd,'rotated.mat'));
%             disp('Rotated data already found in rotated.mat')
%             % Run the rotation anyway if manual override is specified
%             if manual_rot_overwrite == 1
%                 [rot_x,rot_y,rotang] = rotate_traj(x,y);
%             end
%         catch
%             [rot_x,rot_y,rotang] = rotate_traj(x,y);
%         end
%     end
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
    goal_R = 1/4*range(x) + min(x);

    %% Establish maze arm boundaries.
%     w = (ymax-ymin)/7.9; %in cm;   Width of arms. Currently getting wider than actual?
%     l = (xmax-xmin)/1.3; %80; Length of center stem (minus rest and end areas).
%     rest = (xmax-xmin)/8; % rest box
%     %Find center arm borders.
%     if exist('centerarm_manual.mat', 'file')
%         load('centerarm_manual.mat', 'center', 'xmin', 'xmax', 'ymin', 'ymax');
%         disp('Loading manually entered center arm location')
%         w = (ymax-ymin)/5; %40;   Width of arms.
%         l = ((xmax-xmin)-(center.x(2) - center.x(1)))/2;
%     else
%        center = getcenterarm_RW(x,y,w,l,rest);
%    end

% lower boundary, upper boundary for y
% left boundary, right boundary for x

    center.y = [rest_y,choice_enter];
    center.x = [center_L,center_R]; % why are they repeated again?
    
    % left arm (between choice point and goal)
    approach_l.x = [goal_L,choice_exitL];
    approach_l.y = [ymax,choice_enter];

    % right arm (between choice point and goal)
    approach_r.x = [goal_R,choice_exitR,choice_exitR,goal_R];
    approach_r.y = [ymax,ymax,choice_enter,choice_enter];

    %left return
    return_l.x = [center_L,xmin,xmin,center_L];
    return_l.y = [rest_y,rest_y,choice_enter,choice_enter];

    %Right return
    return_r.x = [center_R,xmax,xmax,center_R];
    return_r.y = [rest_y,rest_y,choice_enter,choice_enter];
  
    %left goal
    goal_l.x = [goal_L,xmin,xmin,goal_L];
    goal_l.y = [choice_enter,choice_enter,ymax,ymax];

    %right goal
    goal_r.x = [goal_R,xmax,xmax,goal_R];
    goal_r.y = [choice_enter,choice_enter,ymax,ymax];
    
    %Choice
    choice.x = [choice_exitL,choice_exitR,choice_exitR,choice_exitL];
    choice.y = [choice_enter,ymax,ymax,choice_enter];

%     %Left arm.
%     left.x = [xmin+l, xmax, xmax, xmin+l];
%     left.y = [ymin, ymin, ymin+w, ymin+w];     
%     %Right arm.
%     right.x = left.x;
%     right.y = [ymax-w, ymax-w, ymax, ymax];   
%     %Left return.
%     return_l.x = [xmax-l, xmax, xmax, xmax-l];
%     return_l.y = [ymin+w, ymin+w, center.y(1), center.y(1)];     
%     %Right return.
%     return_r.x = return_l.x;
%     return_r.y = [center.y(3), center.y(3), ymax-w, ymax-w];    
%     %Choice.
%     choice.x = [xmin, xmin+l, xmin+l, xmin];
%     choice.y = [center.y(1), center.y(1), center.y(3), center.y(3)];     
%     %Left approach.
%     approach_l.x = choice.x;
%     approach_l.y = [ymin, ymin, center.y(1), center.y(1)];     
%     %Right approach.
%     approach_r.x = choice.x;
%     approach_r.y = [center.y(3), center.y(3), ymax ymax];    
%     %Base.
%     base.x = return_l.x;
%     base.y = choice.y;     
%     %Right Goal
%     xmin_g = 0.9*(xmax-xmin)+xmin; %0.7 before
%     xmax_g = xmin_g + 0.1*(xmax-xmin);
%     goal_r.x = [xmin_g xmax_g xmax_g xmin_g]; % Seems to work ok for our current maze...
%     goal_r.y = right.y;     
%     %Left Goal
%     goal_l.x = goal_r.x;
%     goal_l.y = left.y;
    
    %% Check with plot.
    if plot_check == 1 % Skip plotting if skip_rot_check = 1;
        figure(555);
        plot(x,y);
        hold on;
        plot([left.x left.x(1)],[left.y left.y(1)], 'k-', ...
            [right.x right.x(1)],[right.y right.y(1)], 'k-', ...
            [return_l.x return_l.x(1)], [return_l.y return_l.y(1)], 'k-', ...
            [return_r.x return_r.x(1)], [return_r.y return_r.y(1)], 'k-', ...
            [center.x center.x(1)], [center.y center.y(1)], 'k-', ...
            [base.x base.x(1)], [base.y base.y(1)], 'g-', ...
            [approach_l.x approach_l.x(1)], [approach_l.y approach_l.y(1)], 'k-', ...
            [approach_r.x approach_r.x(1)], [approach_r.y approach_r.y(1)], 'k-', ...
            [choice.x choice.x(1)], [choice.y choice.y(1)], 'r-');
        hold off;
    end
%         %Sanity check for trajectory rotation.
%         if manual_rot_overwrite == 1
%             satisfied = input('Are you satisfied with the rotation? Enter y or n-->','s');
%         elseif manual_rot_overwrite == 0
%             satisfied = 'y';
%         end
%         
%         if strcmp(satisfied,'y')       %Break.
%             skewed = 0;
%             if manual_rot_overwrite == 1
%                 save rotated rotang rot_x rot_y;
%             end
%         elseif strcmp(satisfied,'n')  %Delete last rotation and try again.
%             if exist(fullfile(pwd, 'rotated.mat'), 'file') == 2
%                 delete rotated.mat;
%             end
%             close all;
%         end
%     elseif skip_rot_check == 1
%         skewed = 0;
%     end
%     
%end

%% Output. 
    %bounds.base = base; 
    bounds.center = center; 
    bounds.choice = choice;
    bounds.approach_l = approach_l;
    bounds.approach_r = approach_r; 
    %bounds.left = left; 
    %bounds.right = right; 
    bounds.return_l = return_l;
    bounds.return_r = return_r; 
    bounds.goal_l = goal_l;
    bounds.goal_r = goal_r;
end
