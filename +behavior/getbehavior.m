function optitrack = getbehavior(varargin)
% firstimport the raw behavioral tracking into a Matlab struct:
% the optitrack output contains timestamps (from optitrack) and position data 
%    optitrack.timestamps      timestamps from optitrack
%    optitrack.position.x      x-position
%    optitrack.position.y      y-position
%    optitrack.position.z      z-position
%    optitrack.speed           speed
%    optitrack.sr              samplingrate

% then add interpolated data points, as well as rotate if necessary.

%% Get inputs
p = inputParser;
addRequired(p,'basepath',@ischar);
addRequired(p,'basename',@ischar);
addRequired(p,'behaviorpath',@ischar);
addParameter(p,'saveMat',true,@islogical);
addParameter(p,'rotateposition',true,@islogical);
addParameter(p,'degreesrotation',180,@isnumeric);
addParameter(p,'centeronzero',true,@islogical);

parse(p,varargin{:});
basepath = p.Results.basepath;
basename = p.Results.basename;
behaviorpath = p.Results.behaviorpath;
saveMat = p.Results.saveMat;
rotateposition = p.Results.rotateposition;
degreesrotation = p.Results.degreesrotation;
centeronzero = p.Results.centeronzero;

%% Get filenames

behdir = dir([behaviorpath '*.csv']);
for f = 1:length(behdir)
    optifilenames{f} = ['behaviorfiles/' behdir(f).name];
end

%RW version has correct cell references, making sure speed/acc aren't nans,
%and then fixing the cm conversion (keeping as is, in cm, rather than *100)
optitrack = optitrack2buzcode_RW('basepath',basepath,'basename',basename,...
    'filenames',optifilenames,saveMat=false);

%don't use the internal saveMat cause the session thing throws it
% 
% 
%     hBrushing    = findall(gca,'tag','Brushing');
%     brushData    = get(hBrushing, {'Xdata','Ydata'});
%     brushIdx     = ~isnan(brushData{1}); % Puts a '1' at indices that aren't NAN
%     brushXData   = brushData{1}(brushIdx);
%     brushYData   = brushData{2}(brushIdx);
%     
%% Interpolate
% Interp
rawx = optitrack.position.x;
rawy = optitrack.position.y;
rawz = optitrack.position.z;

x1 = 1:numel(rawx);
%// Indices of NaNs
t2 = find(~isnan(rawx));
%// Replace NaNs with the closest non-NaNs
optitrack.position.interpolatedx = interp1(x1(t2),rawx(t2),x1,'nearest'); % interpolated values!
optitrack.position.interpolatedy = interp1(x1(t2),rawy(t2),x1,'nearest'); % interpolated values
optitrack.position.interpolatedz = interp1(x1(t2),rawz(t2),x1,'nearest');

optitrack.position.interpolatedx(1)=optitrack.position.interpolatedx(2);
optitrack.position.interpolatedy(1)=optitrack.position.interpolatedy(2);
optitrack.position.interpolatedz(1)=optitrack.position.interpolatedz(2);

if rotateposition == 1
    x = optitrack.position.interpolatedx;
    y = optitrack.position.interpolatedy;
    optitrack.position.interpolatedx = x*cosd(degreesrotation) - y*sind(degreesrotation);
    optitrack.position.interpolatedy = x*sind(degreesrotation) + y*cosd(degreesrotation);
end

if centeronzero == 1
    %maze.position.interpolatedx = maze.position.interpolatedx - min(maze.position.interpolatedx);
    %maze.position.interpolatedy = maze.position.interpolatedy - min(maze.position.interpolatedy);

    xmean = mean(optitrack.position.interpolatedx,'omitnan');
    ymean = mean(optitrack.position.interpolatedy,'omitnan');
    optitrack.position.interpolatedx = optitrack.position.interpolatedx - xmean;
    optitrack.position.interpolatedy = optitrack.position.interpolatedy - ymean;

end

if saveMat == true
    save(fullfile(basepath,[basename '.dat.optitrack.behavior.mat']),'optitrack')
end

