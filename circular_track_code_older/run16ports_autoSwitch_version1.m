function [] = run16ports_autoSwitch_version1()

% detect which sensors have been activated and deliver water in response
% also calculates performance of rat, as well as specific inside/outside
% performance (specific to rachel black/white circular maze task)
% prefill and fills opposite well rather than current well nais -BK
% testwells(0,.75)
% testwells(1,.75)
% testwells(2,.75)
% testwells(3,.75)

%%  Target pairs + general params

% rat = 'test';

% sensor names from the PCB BOARD/PHYSICAL SENSORS, [# #], double
% sensors_pcb{2} = [12 13]; % pair 2
% valves_pcb{2} = [1 2];
 % sensors_pcb{1} = [2 3]; % pair 1
 % valves_pcb{1} = [3 0];1,
% Brian stole some of your code and now this is for him % i'm furious how dare you
rat = 'BK_Creampuff_watertraining_20240624';   % Wow brian stole some of your code here too teehee % this is the LAST STRAW
sensors_pcb{1} = [2 3]; 
valves_pcb{1} = [0 2];

%thresh = .85; % 85% performance
pumpOpen = .1;% .04; % in seconds, length of time water is delivered (small because of the gravity system)
%don't go lower than 0.02! solenoid won't work. sometimes unreliable up to 0.07 as well.

%% don't edit these param values:

bInside_pcb = [12 13 14 15];
bOutside_pcb = [4 5 6 7];
wInside_pcb = [8 9 10 11];
wOutside_pcb = [0 1 2 3];

bInside = convertSensorName(bInside_pcb); % spits out CELL index
bOutside = convertSensorName(bOutside_pcb); % spits out CELL index
wInside = convertSensorName(wInside_pcb); % spits out CELL index
wOutside = convertSensorName(wOutside_pcb); % spits out CELL index

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

%% initialize spacebar detection

spacefig = figure('Name','I am just here to detect space bar - keep me in focus');
set(spacefig,'KeyPressFcn',@isspacebar); %see below')
drawnow;

%% initialize file

starttime = datetime;
%lastPrintTime = starttime;
t = replace(char(timeofday(starttime)),':','.');
%ztreference = datetime([date '12:00:00:000'], 'Format','dd-mmm-yyyy hh:mm:ss:SSS')
filename = ['E:\Rachel\Datapoints\datapoints_' num2str(yyyymmdd(datetime)) '_'  t '_' rat '.txt'];
fID = fopen(filename,'a');

%% Create cleanup obj to save figs + close file upon ctrl+c
di = [];
do = [];
cleanup = onCleanup(@()myCleanupFun(h,fID,rat,di,do));

%% for pair 1 and 2

for pair = 1:2

    disp(['pair = ' num2str(pair)])
    clear di do
    %% receiving input from national instruments
    % To read inputs from 16 sensors on NI board and output to up to 4 solenoids

    di = DigitalInput(); %from 16 sensors
    do = DigitalOutput(); %out to up to 4 solenoids

    % Convert sensor names. send in as double, comes out as cell w/ NI name
    sensors = convertSensorName(sensors_pcb{pair}); % spits out CELL index
    valves = convertSensorName(valves_pcb{pair}); % spits out CELL index

    %% prefill wells
      testwells(valves_pcb{pair}(1),pumpOpen);
      testwells(valves_pcb{pair}(2),pumpOpen);

    %% RUN

    lastSensor = 0;
    lastCorrect = 0;
    performanceVec = [];
    bwVec = [];
    n = 1;
    newsensortime = [];
    tri = 0;
    switchpairs = 0;

    while true && switchpairs == 0

        drawnow;

        %%% switch condition %%%
        % read isspacebar below! switches to pair 2 with spacebar push
        %%%%%%%%%%%%%%%%%%%%%%%%

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
            do.toggleNI(valves(2),pumpOpen) % output,time(seconds)
            lastCorrect = sensors(1);
            firstport = sensors(1);
            outcome = 1;
            %%%% triggered and correct
        elseif activeport == sensors(1) && lastCorrect == sensors(2)
            do.toggleNI(valves(2),pumpOpen) % output,time(seconds)
            lastCorrect = sensors(1);
            outcome = 1;
            %%%% triggered but wrong
        elseif activeport == sensors(1) && lastCorrect ~= 0 && lastCorrect ~= sensors(2)
            outcome = 0;

            %% Port 2 triggered
            %%%% first port in session
        elseif activeport == sensors(2) && lastCorrect == 0
            do.toggleNI(valves(1),pumpOpen) % output,time(seconds)
            lastCorrect = sensors(2);
            firstport = sensors(2);
            outcome = 1;

            %%%% triggered and correct
        elseif activeport == sensors(2) && lastCorrect == sensors(1)
            do.toggleNI(valves(1),pumpOpen) % output,time(seconds)
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
        if repeat == 0 || seconds(d-newsensortime)>= 15 % %if at least 15 seconds since the rat initiated new sensor visit
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

            drawnow; % to force fig update

            n = n+1;
            %% Write to file and send to open ephys

            %if seconds(d-lastPrintTime) >= 5 % if at least 5 seconds since the last print time
            fprintf(fID,' %s %d %i \n', string(d,'dd-MMM-yyyy hh:mm:ss:SSS'), activeport, performance); %then write to file
            lastPrintTime = datetime; lastPrintTime.Format = 'hh:mm:ss:SSS';

            % send to open ephys
            %end

            if outcome == 1 && activeport ~= firstport % just to print # of trials
                tri = tri + 1;
                disp(['Trial ' num2str(tri) ' complete'])
            end

        end %  if repeat == 0 || seconds(d-newsensortime)>= 15
    end % while true
end % for pair = 1:2
end %run16ports_autoSwitch

function myCleanupFun(h,fID,rat,di,do)
% so when ctrl+c it closes the file and saves the figures
fclose(fID);
d = datetime;
saveas(h, ['E:\Brian\WaterTrainingFigures\' char(d,'dd-MMM-yyyy_hh.mm.ss.SSS') '_' rat '_performance.fig'])
saveas(h, ['E:\Brian\WaterTrainingFigures\' char(d,'dd-MMM-yyyy_hh.mm.ss.SSS') '_' rat '_performance.pdf'])
delete(di)
delete(do)
disp('File and figures saved') %if you close fig before saving you'll get saveas error
end

function isspacebar(src,event)
if strcmp(event.Key,'space')
    assignin("caller","switchpairs",1)
 %   xline()
  %  switchpairs = 1; %switch to pair 2
    disp('Switching to pair 2: RING BELL')
end
end
