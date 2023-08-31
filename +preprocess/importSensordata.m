function sensordata = importSensordata(foldername)
% To import sensor data from matlab
% foldername = '/data/ExperimentsRSW/CircularMaze/20230614/merged_20230614/sensors/'

files = dir([foldername '*.txt']); % get all sensor files

for f = 1:length(files)
    thefile = fullfile(foldername,files(f).name);
    fID = fopen(thefile);
    file = fread(fID,[3,Inf]);








































end