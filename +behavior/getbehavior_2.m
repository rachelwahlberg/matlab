function optitrack = getbehavior(varargin)

%% Get inputs
p = inputParser;
addRequired(p,'basename',@ischar);
addRequired(p,'basepath',@ischar);
addRequired(p,'behaviorpath',@ischar);
addParameter(p,'saveMat',true,@islogical);
addParameter(p,'rotateposition',true,@islogical);
addParameter(p,'degreesrotation',-90,@isnumeric);

parse(p,varargin{:});
basename = p.Results.basename;
basepath = p.Results.basepath;
behaviorpath = p.Results.behaviorpath;
saveMat = p.Results.saveMat;
rotateposition = p.Results.rotateposition;
degreesrotation = p.Results.degreesrotation;

%% Get filenames

behdir = dir([behaviorpath '*.csv']);
for f = 1:length(behdir)
    optifilenames{f} = ['behaviorfiles/' behdir(f).name];
end
optitrack = optitrack2buzcode_RW('basepath',basepath,'basename',basename,...
    'filenames',optifilenames,saveMat=false);
%don't use the internal saveMat cause the session thing throws it

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
    optitrack.position.interpolatedx = x*cosd(optitrack.rotation) - y*sind(optitrack.rotation);
    optitrack.position.interpolatedy = x*sind(optitrack.rotation) + y*cosd(optitrack.rotation);
end

if saveMat == true
    save(fullfile(basepath,[basename '.dat.optitrack.behavior.mat']),'optitrack')
end

