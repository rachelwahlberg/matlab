function optitrack = outofmaze(varargin)
% to skip over this function set 'checkforoutofmaze' as false

p = inputParser;
addRequired(p,'optitrack',@isstruct);
addParameter(p,'checkforoutofmaze',true,@islogical)
addParameter(p,'basepath',[],@ischar)
addParameter(p,'basename',[],@ischar)

parse(p,varargin{:})
optitrack = p.Results.optitrack;
checkforoutofmaze = p.Results.checkforoutofmaze;
basepath = p.Results.basepath;
basename = p.Results.basename;

%to use:
% will plot a figure, if you want to remove any points type 1, then enter
% will be paused. use the brushing tool to select points to remove, and
% then set those values to nan (in the figure gui under tools)
% press enter to have those nan values incorporated into optitrack struct.

if checkforoutofmaze == 1
    figure
    h = scatter(optitrack.position.interpolatedx,optitrack.position.interpolatedy);

    prompt = "Do you want to remove out of maze points? 1 or 0";
    remove = logical(input(prompt)); %
    
    if remove
        disp(['use the brushing tool to select points to remove, and' ...
            'then set those values to nan (in the figure gui under tools)' ...
            'press enter to have those nan values incorporated into optitrack struct.'])
        pause
        % use the brush tool to select points to remove, and then set these values
        % to nans.
        for s = 1:length(h.YData)
            if isnan(h.XData(s))
                h.YData(s) = nan;
            end

            if isnan(h.YData(s))
                h.XData(s) = nan;
            end
        end

        optitrack.position.interpolatedx = h.XData;
        optitrack.position.interpolatedy = h.YData;

        if ~isempty(basepath) && ~isempty(basename)
        save(fullfile(basepath, [basename '.dat.optitrack.behavior.mat']),'optitrack');
        end

        figure
        hold all
        title('adjusted figure')
        plot(optitrack.position.interpolatedx,optitrack.position.interpolatedy)

    else
        disp("Not removing any out of maze points")
    end
end

