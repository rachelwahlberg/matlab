function [] = sensorData(varargin)
% detect which sensors have been activated and deliver water in response
% also calculates performance of rat, as well as specific inside/outside
% performance (specific to rachel black/white circular maze task)

%% Initial params

p = inputParser;
addRequired(p,'di',@isobject)
addRequired(p,'do',@isobject)
addRequired(p,'sensors',@isnumeric)
addRequired(p,'valves',@isnumeric)
addParameter(p,'openTime',.15,@isnumeric) % valve open time in s
addParameter(p,'calculateInsideOutside',false,@islogical)
addParameter(p,'bInside',[],@isnumeric)
addParameter(p,'wInside',[],@isnumeric)
addParameter(p,'bOutside',[],@isnumeric)
addParameter(p,'wOutside',[],@isnumeric)

parse(p,varargin{:});
di = p.Results.di;
do = p.Results.do;
sensors = p.Results.sensors;
valves = p.Results.valves;
openTime = p.Results.openTime;
calculateInsideOutside = p.Results.calculateInsideOutside;
bInside = p.Results.bInside;
wInside = p.Results.wInside;
bOutside = p.Results.bOutside;
wOutside = p.Results.wOutside;

starttime = datetime;
t = replace(char(timeofday(starttime)),':','.');
filename = ['datapoints_' num2str(yyyymmdd(datetime)) '_' t '.txt'];
fID = fopen(filename,'a');

lastSensor = 0;
totalcorrect = 0;
bTotal = 0;
wTotal = 0;
n = 1; bN = 1; wN = 1;

%% initialize figure

h = figure;
hold all
h1 = subplot(1,2,1,'Parent',h);
hold on
title('Overall Performance') % how often they go to the specific required port
xlabel('Minutes from start')
ylabel('Performance')
h2 = subplot(1,2,2,'Parent',h);
hold on
title('Inside/Outside Performance') % how often they sample a port on the correct side (inside on b, outside on w)
xlabel('Minutes from start')
ylabel('Performance')

%% Run while loop 

while true

    port =di.readNI; % default is one if it's not being triggered, zero if it is.
    activeport=find(port==1);
    if isempty(activeport)
        activeport = 0;
        continue
    end

    %% Port 1 triggered

    %%%% first port in session
    if activeport == sensors(1) && lastSensor == 0
        do.toggle(valves(1),openTime) % output,time(seconds) 
        lastSensor = 1;
        outcome = 1;
    end

    %%%% triggered and correct
    if activeport == sensors(1) && lastSensor == sensors(2)
        do.toggle(valves(1),openTime) % output,time(seconds) 
        lastSensor = 1;
        outcome = 1;
    end

    %%%% triggered but wrong
    if activeport == sensors(1) && lastSensor ~= 0 && lastSensor ~= sensors(2)
        outcome = 0;
    end

    %% Port 2 triggered

    %%%% first port in session
    if activeport == sensors(2) && lastSensor == 0
        do.toggle(valves(2),openTime) % output,time(seconds) 
        lastSensor = 2;
        outcome = 1;
    end

    %%%% triggered and correct
    if activeport == sensors(2) && lastSensor == sensors(1)
        do.toggle(valves(2),openTime) % output,time(seconds) 
        lastSensor = 2;
        outcome = 1;
    end

    %%%% triggered but wrong
    if activeport == sensors(2) && lastSensor ~= 0 && lastSensor ~= sensors(1)
        outcome = 0;
    end

    %% Any other ports triggered

    if (activeport ~= sensors(1) || activeport ~= sensors(2)) ...
            && activeport ~= 0
        outcome = 0;
    end

    %% Calculate performance
    totalcorrect = totalcorrect + outcome;
    performance = totalcorrect/n;
    n = n+1;
    d  = datetime;
    scatter(h1,minutes(d-starttime),performance)
    disp(['Performance = ' num2str(performance)])

    %% Calculate b/w performance 
    if calculateInsideOutside
        if any(bInside) == activePort
            bOutcome = 1;
        end

        if any(bOutside) == activePort
            bOutcome = 0;
        end

        if any(wInside) == activePort
            wOutcome = 0;
        end

        if any(wOutside) == activePort
            wOutcome = 1;
        end

       if any([bInside bOutside] == activePort)
           bTotal = bTotal + bOutcome;
           bPerformance = bTotal/bN;
           bN = bN + 1;
           scatter(h2,minutes(d-starttime),bPerformance,'black','filled')
       elseif any([wInside wOutside] == activePort)
           wTotal = wTotal + wOutcome;
           wPerformance = wTotal/wN;
           wN = wN + 1;
           scatter(h2,minutes(d-starttime),wPerformance,'black')
       end
    end


    %% Write to file and send to open ephys

     fprintf(fID,'%d \n', activeport);

     % send to open ephys
end




end