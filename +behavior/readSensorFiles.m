function [day,time,port] = readSensorFiles(filepath)

%filepath = '/data/20230614/sensorData/';

day = [];
time = [];
port = [];

allfiles = dir([filepath '*.txt']);

for f = 1:length(allfiles)
    fID = fopen([filepath allfiles(f).name],'r');
    file = textscan(fID,'%s %s %d');
    fclose(fID);

    if isempty(day)
        day = file{1}{1};
    end

    time = [time; file{2}];
    port = [port; file{3}];
end

t = datetime(time{1},InputFormat="hh:mm:ss:SSS",'hh:mm:ss:SSS')







