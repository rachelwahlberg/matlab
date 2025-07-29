function [] = run16ports()

% detect which sensors have been activated and deliver water in response
% also calculates performance of rat, as well as specific inside/outside
% performance (specific to rachel black/white circular maze task)

%%  Target pairs + general params

rat = 'Harry_DG_session1';
%thresh = .85; % 85% performance
pumpOpen = .01; % in seconds, length of time water is delivered (small because of the gravity system)

% sensor names from the PCB BOARD/PHYSICAL SENSORS, [# #], double
%sensors_pcb = [0 13];
%valves_pcb = [1 2];
sensors_pcb = [2 14];
valves_pcb = [0 3];
bInside_pcb = [12 13 14 15];
bOutside_pcb = [4 5 6 7];
wInside_pcb = [8 9 10 11];
wOutside_pcb = [0 1 2 3];

% Convert sensor names. send in as double, comes out as cell w/ NI name
sensors = convertSensorName(sensors_pcb); % spits out CELL index
valves = convertSensorName(valves_pcb); % spits out CELL index
bInside = convertSensorName(bInside_pcb); % spits out CELL index
bOutside = convertSensorName(bOutside_pcb); % spits out CELL index
wInside = convertSensorName(wInside_pcb); % spits out CELL index
wOutside = convertSensorName(wOutside_pcb); % spits out CELL index

%% receiving input from national instruments
% To read inputs from 16 sensors on NI board and output to up to 4 solenoids

di = DigitalInput();
do = DigitalOutput();

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

%% initialize file

starttime = datetime;
lastPrintTime = starttime;
t = replace(char(timeofday(starttime)),':','.');
filename = ['Datapoints/datapoints_' num2str(yyyymmdd(datetime)) '_'  t '_' rat '_.txt'];
fID = fopen(filename,'a');

%% Create cleanup obj to save figs + close file upon ctrl+c

cleanup = onCleanup(@()myCleanupFun(h,fID,rat));

%% RUN

lastSensor = 0;
lastCorrect = 0;
performanceVec = [];
bwVec = [];
n = 1;
newsensortime = [];
% 
% for i=1:10
%   do.toggleNI(valves(1),pumpOpen) % output,time(seconds)
%   pause(2)
% end
% 
while true

    port =di.readNI; % default is zero if it's not being triggered, one if it is.
    activeport=find(port==1);
 %   activeport

    if isempty(activeport)
        activeport = 0;
        continue
    end

    %% Port 1 triggered

    %%%% first port in session
    %identify if it's a repeat port
    if activeport == lastSensor
        repeat = 1;
    else
        repeat = 0;
        newsensortime = datetime;
    end

    %%%% water trigger rules
    if activeport == sensors(1) && lastCorrect == 0
        do.toggleNI(valves(1),pumpOpen) % output,time(seconds)
       lastCorrect = sensors(1);
        outcome = 1;
        %%%% triggered and correct
    elseif activeport == sensors(1) && lastCorrect == sensors(2)
        do.toggleNI(valves(1),pumpOpen) % output,time(seconds)
       lastCorrect = sensors(1);
        outcome = 1;
        %%%% triggered but wrong
    elseif activeport == sensors(1) && lastCorrect ~= 0 && lastCorrect ~= sensors(2)
        outcome = 0;

        %% Port 2 triggered
        %%%% first port in session
    elseif activeport == sensors(2) && lastCorrect == 0
        do.toggleNI(valves(2),pumpOpen) % output,time(seconds)
       lastCorrect = sensors(2);
        outcome = 1;

        %%%% triggered and correct
    elseif activeport == sensors(2) && lastCorrect == sensors(1)
        do.toggleNI(valves(2),pumpOpen) % output,time(seconds)
        lastCorrect = sensors(2);
        outcome = 1;

        %%%% triggered but wrong
    elseif activeport == sensors(2) && lastCorrect ~= 0 && lastCorrect ~= sensors(1)
        outcome = 0;

        %% Any other ports triggered
    
    elseif activeport ~= sensors(1) && activeport ~= sensors(2) ...
            && activeport ~= 0
        outcome = 0;
    end

    lastSensor = activeport;

    %% Calculate general performance
    d= datetime;
    if repeat == 0 || seconds(d-newsensortime)>= 15 %if at least 15 seconds since the rat initiated new sensor visit
        if repeat == 1
            newsensortime = datetime;
        end

        performanceVec(n) = outcome;
        if n <10
            performance = sum(performanceVec)/length(performanceVec);
        else
            performance = sum(performanceVec(n-9:n))/10;
        end
        scatter(h1,minutes(d-starttime),performance)
        disp(['Performance = ' num2str(performance)])

        %% Calculate b/w performance

        if any(bInside == activeport)
            bwVec(n) = 1;
        elseif any(bOutside == activeport)
            bwVec(n) = 2;
        elseif any(wInside == activeport)
            bwVec(n) = 3;
        elseif any(wOutside == activeport)
            bwVec(n) = 4;
        end
        
        if any([bInside bOutside] == activeport)
            if n < 10
                bPerformance = sum(bwVec == 1)/(sum(bwVec == 1) + sum(bwVec == 2));
            else
                bPerformance = sum(bwVec(n-9:n) == 1)/(sum(bwVec(n-9:n) == 1) + sum(bwVec(n-9:n) == 2));
            end

            scatter(h2,minutes(d-starttime),bPerformance,'black','filled')
            disp(['bPerformance = ' num2str(bPerformance)])
        elseif any([wInside wOutside] == activeport)
            
            if n < 10
                wPerformance = sum(bwVec == 4)/(sum(bwVec == 3) + sum(bwVec == 4));
            else
                wPerformance = sum(bwVec(n-9:n) == 4)/(sum(bwVec(n-9:n) == 3) + sum(bwVec(n-9:n) == 4));
            end

            scatter(h2,minutes(d-starttime),wPerformance,'black')
            disp(['wPerformance = ' num2str(wPerformance)])
        end

        n = n+1;
        %% Write to file and send to open ephys

        if seconds(d-lastPrintTime) >= 5 % if at least 5 seconds since the last print time
            fprintf(fID,' %s %d \n', string(d,'dd-MMM-yyyy hh:mm:ss:SSS'), activeport); %then write to file
            lastPrintTime = datetime; lastPrintTime.Format = 'hh:mm:ss:SSS';

            % send to open ephys
        end

    end
end

end

    function myCleanupFun(h,fID,rat)
        % so when ctrl+c it closes the file and saves the figures
        fclose(fID);
        d = datetime;
        saveas(h, ['Figures/' char(d,'dd-MMM-yyyy_hh.mm.ss.SSS') '_' rat '_performance.fig'])
        saveas(h, ['Figures/' char(d,'dd-MMM-yyyy_hh.mm.ss.SSS') '_' rat '_performance.pdf'])
        disp('File and figures saved') %if you close fig before saving you'll get saveas error
    end
