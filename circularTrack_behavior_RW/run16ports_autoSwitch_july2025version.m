function [] = run16ports_autoSwitch_july2025version()

% for AUTOMATIC water delivery. See run16ports_autoSwitch for auto switching using sensors.
% detect which sensors have been activated and deliver water in response
% also calculates performance of rat, as well as specific inside/outside
% % performance (specific to rachel black/white circular maze task)

% testwells(0,0.5)
% testwells(1,0.5)
% testwells(2,0.5)
% testwells(3,0.5)

 %%% IMPORTANT:
 % to close and save, hit return, and then ctrl+c or quit. HItting return
 % is VERY IMPORTANT because it initalizes the cleanup func so that all the
 % data gets saved (see the keyboard control function). Can't initialize outside of the while loop because the cleanup func uses
 % the variables defined at the time of initalization.

 % ALSO if you say "session 1" twice in a day or something the data will be
 % appended to the original session 1 file and will overwrite the old fig. so be careful

 % IF RECORDING, run this before your pre session to initialize the
 % directories before setting your ephys/motive folders.

%% Initialize session and/or rat folder/csv file
% to run a test name the rat Test and session type Test. Will go into a
% separate folder to keep actual folders clean.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% RUN THIS SECTION BEFORE SETTING EPHYS/MOTIVE PATHS AT BEGINNING OF RECORDING DAY. %%%% 

rats_basepath = 'E:\Rachel\CircularTrack\Recording_Rats'; % Rat will be appended on
currentdate =  char(datetime('today','format','yyyy-MM-dd'));

%input session info
prompt = {'Rat name (uppercase) ("Test" for dryrun):','Session Type (WaterTraining, TaskTraining, Recording, or TEST)', 'Session number:', 'Region (CA1 or DG):'};
dlgtitle = 'Session information';
fieldsize = [1 45; 1 45; 1 45; 1 45]; 
answer = inputdlg(prompt,dlgtitle,fieldsize);

rat = answer{1};
session_type = answer{2};
session_number = answer{3}; %within a day, so this'll really only be 1 or 2
region = answer{4};

if strcmp(session_type,'WaterTraining') ~= 1 && strcmp(session_type,'TaskTraining') ~= 1 && strcmp(session_type,'Recording') ~= 1 && strcmp(session_type,'TEST') ~= 1
    disp('Session type must be WaterTraining, TaskTraining, Recording, or TEST');
    return
end

% creates folders for today's date if they don't already exist
session_basepath = initializeDirectories(rats_basepath,rat,session_type,region,currentdate);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% don't edit these param values:

bInside_pcb = [12 13 14 15];
bOutside_pcb = [4 5 6 7];
wInside_pcb = [8 9 10 11];
wOutside_pcb = [0 1 2 3];

bInside = convertSensorName(bInside_pcb); % spits out CELL index
bOutside = convertSensorName(bOutside_pcb); % spits out CELL index
wInside = convertSensorName(wInside_pcb); % spits out CELL index
wOutside = convertSensorName(wOutside_pcb); % spits out CELL index

%sound file 
bellSoundPath = 'sounds/service-bell_daniel_simion.mp3'; %in PhoPlaySoundsMatlab

%% session checks
%because I have a flighty memory okay

if strcmp(answer{2},'Recording') == 1
    prompt = {'Did you set your Motive path? (y/n)','Did you set your OpenEphys path? (y/n)', 'Did you start the video? (y/n)', ...
        'Did you alcohol swab the ports? (y/n)','Did you start the recording? (=^_^=) (y/n)'};
    dlgtitle = 'tis for your own good';
    fieldsize = [1 60; 1 60; 1 60; 1 60; 1 60;];
    answer  = inputdlg(prompt,dlgtitle,fieldsize);

    if strcmp(answer{1},'y') ~= 1 || strcmp(answer{2},'y') ~= 1 || strcmp(answer{3},'y') ~= 1 || strcmp(answer{4},'y') ~= 1 || strcmp(answer{5},'y') ~= 1
        disp('come on, son')
        return
    end
end

%%  Target pairs + general params
% sensor names from the PCB BOARD/PHYSICAL SENSORS, [# #], double
sensors_pcb{1} = [12 13]; % 12 13
valves_pcb{1} = [2 1]; 

sensors_pcb{2} = [14 15]; %14 15
valves_pcb{2} = [0 3];

pumpOpen = 0.08; % .04; % in seconds, length of time water is delivered (small because of the gravity system)
%don't go lower than 0.02! solenoid won't work. sometimes unreliable up to 0.07 as well.

% bInside_pcb = [12 13 14 15];
% bOutside_pcb = [4 5 6 7];
% wInside_pcb = [8 9 10 11];
% wOutside_pcb = [0 1 2 3];
% 
% %% initialize figure
% figure
% h = gca;
% hold all
% title('Water delivery times')
% xlabel('Minutes from start')

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

spacefig = figure('Name','I am just here to detect keyboard - keep me in focus');
set(spacefig,'KeyPressFcn',@keyboardcontrol); %see below')
drawnow;

%% initialize files

starttime = datetime;
starttime.Format = 'HH:mm:ss:SSS';
disp(strcat(rat," file start time is ",string(starttime)))

%%%%%% port visits %%%%%%
fnamePorts = [session_basepath '\National_Instruments\portTimestamps_auto_' rat '_' region '_' currentdate '_session' session_number '.txt'];
fIDd = fopen(fnamePorts,'a');

%%%%% start/transition times %%%%%
fnameTransition = [session_basepath '\National_Instruments\portTransitions_auto_' rat '_' region '_' currentdate '_session' session_number '.txt'];
fIDt = fopen(fnameTransition,'a');
fprintf(fIDt,' %s %i %i \n', string(starttime),sensors_pcb{1}(1),sensors_pcb{1}(2)); %then write to file

% %%%%% performance %%%%%
% fnamePerformance = [session_basepath '\National_Instruments\performance_auto_' rat '_' region '_' currentdate '_session' session_number '.txt'];
% fIDp = fopen(fnamePerformance,'a');
% %fprintf(fIDp,' %s %i %i \n', string(starttime),sensors_pcb{1}(1),sensors_pcb{1}(2)); %then write to file

%%%% consolidated csv over sessions, w/ transition times and notes %%%%%
%will be created or appended in the cleanup func

headers = {'Date','Session','Session_Type','Time_Started','Time_Switched','Time_Ended','Ports_1' 'Ports_2','NumberTrials_1','NumberTrials_2','Notes'};
sessionTable = cell2table(cell(1,11),'VariableNames',headers); %initialize
sessionTable.Date = datetime('today','format','yyyy-MM-dd');
sessionTable.Session = session_number;
sessionTable.Session_Type = session_type;
sessionTable.Time_Started = starttime;
sessionTable.Ports_1 = num2str([sensors_pcb{1}(1),sensors_pcb{1}(2)]);

%% Create cleanup obj to save figs + close file upon ctrl+c, as well as reroute video (if taken) to correct dir
di = [];
do = [];
tri = 0;
closeupshop = 0;
%cleanup = onCleanup(@()myCleanupFun(h,fIDd,fIDt,rat,di,do));

%% for pair 1 and 2

for pair = 1:2
   
    %% receiving input from national instruments
    % To read inputs from 16 sensors on NI board and output to up to 4 solenoids

    disp(['pair = ' num2str(pair)])
    clear do di

    di = DigitalInput(); %from 16 sensors
    do = DigitalOutput(); %out to up to 4 solenoids

    % Convert sensor names. send in as double, comes out as cell w/ NI name
    sensors = convertSensorName(sensors_pcb{pair}); % spits out CELL index
    valves = convertSensorName(valves_pcb{pair}); % spits out CELL index

    if pair == 2 % draw a transition line
        start2 = datetime; start2.Format = 'HH:mm:ss:SSS';
        playLocalAudioFile(bellSoundPath) %make sure speakers are on - indicate that pair has switched with a bell sound.
        disp(strcat(rat," Pair ",num2str(pair)," start time is ",string(start2)))
        xline(h1,minutes(start2-starttime))
         xline(h2  ,minutes(start2-starttime))
        fprintf(fIDt,' %s %i %i \n', string(start2),sensors_pcb{2}(1),sensors_pcb{2}(2)); %then write to file
        
        sessionTable.Time_Switched = start2;
        sessionTable.Ports_2 =  num2str([sensors_pcb{2}(1),sensors_pcb{2}(2)]);
        sessionTable.NumberTrials_1 = tri;
    end

    %% prefill wells
     %testwells(valves_pcb{pair}(1),pumpOpen);
     %testwells(valves_pcb{pair}(2),pumpOpen);

    %% RUN

    lastSensor = 0;
    lastCorrect = 0;
    n = 1;
    performanceVec = [];
    bwVec = []; 
    %newsensortime = [];
    tri = 0;
    switchpairs = 0; 
    activeport = 0;

    while true && switchpairs == 0

        drawnow;

        if closeupshop == 1
             cleanup = onCleanup(@()myCleanupFun(h,fIDd,fIDt,rat,di,do,tri,sessionTable,rats_basepath,session_basepath));
             pause
        end

        %% switch condition %%%
        %read isspacebar below! switches to pair 2 with spacebar push
        %%%%%%%%%%%%%%%%%%%%%%%

        port =di.readNI; % default is zero if it's not being triggered, one if it is.
        activeport=find(port==1);

       % if activeport == 0
       %    continue
       %  end
       if isempty(activeport)
           activeport = 0;
           continue
       end
        
        % get activeport from the keyboard

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
            firstport = sensors(1);
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
            firstport = sensors(2);
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

        %% Write water times to file
        d= datetime;
        if repeat == 0 || seconds(d-newsensortime)>=15
          
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



           %scatter(h,d,1)
            drawnow; % to force fig update
            n = n+1;
            %% Write to file and send to open ephys

            %if seconds(d-lastPrintTime) >= 5 % if at least 5 seconds since the last print time
            fprintf(fIDd,' %s %d %d \n', string(d,'HH:mm:ss:SSS'), activeport-1,performance); %then write to file. -1 to get back to 0 indexing
            lastPrintTime = datetime; lastPrintTime.Format = 'HH:mm:ss:SSS';

            % send to open ephys

            %end

              disp([rat ' visited port ' num2str(activeport - 1) '. Correct: ' num2str(outcome)])

            if outcome == 1 && activeport ~= firstport % just to print # of trials
                tri = tri + 1;
                disp(['Trial ' num2str(tri) ' complete'])
               
            end

            activeport = 0;

        end %  if repeat == 0 || seconds(d-newsensortime)>= 15
    end % while true
end % for pair = 1:2
end %run16ports_autoSwitch

function myCleanupFun(h,fIDd,fIDt,rat,di,do,tri,sessionTable,rats_basepath,session_basepath)
% so when ctrl+c it closes the file and saves the figures

dEnd = datetime;
dEnd.Format = 'HH:mm:ss:SSS';
fprintf(fIDt,' %s \n', string(dEnd)); %write end time to file


prompt = {'Session notes:'};
dlgtitle = 'Session notes';
fieldsize = [20 45]; 
answer = inputdlg(prompt,dlgtitle,fieldsize);

sessionTable.Time_Ended = dEnd;
sessionTable.Notes = answer{1};

if iscell(sessionTable.NumberTrials_1)
    sessionTable.NumberTrials_1 = tri;
else
    sessionTable.NumberTrials_2 = tri;
end


if ~isfile([rats_basepath '\' rat '\' rat '_notes.csv']) %to sythesize across days/sessions
    writetable(sessionTable,[rats_basepath '\' rat '\' rat '_notes.csv'], 'WriteVariableNames',true,'WriteRowNames',true)
else
    writetable(sessionTable,[rats_basepath '\' rat '\' rat '_notes.csv'],'WriteMode','Append','WriteVariableNames',false,'WriteRowNames',true);
end


fclose(fIDd); %datapoints file
fclose(fIDt); %transition file
currentdate =  char(datetime('today','format','yyyy-MM-dd'));
saveas(h, [session_basepath '\Figures\' currentdate '_' rat '_' sessionTable.Session '_waterdeliverytimes.fig'])
%delete(di)
%delete(do)
datetime
disp('File and figures saved') %if you close fig before saving you'll get saveas error
disp('Now go move your video file!')
end

function keyboardcontrol(src,event) %technically manual still works
%MATLAB CODING
if strcmp(event.Key,'space') % switch to pair 2
    assignin("caller","switchpairs",1)
    disp('Switching to pair 2')
% elseif strcmp(event.Key,'q') % visited port 0 
%     assignin("caller","activeport",1)
% elseif strcmp(event.Key,'w') % visited port 1
%     assignin("caller","activeport",2)
% elseif strcmp(event.Key,'e') % visited port 2
%     assignin("caller","activeport",3)
% elseif strcmp(event.Key,'r') % visited port 3
%     assignin("caller","activeport",4)
% elseif strcmp(event.Key,'t') % visited port 4
%     assignin("caller","activeport",5)
% elseif strcmp(event.Key,'y') % visited port 5
%     assignin("caller","activeport",6)
% elseif strcmp(event.Key,'u') % visited port 6
%     assignin("caller","activeport",7)
% elseif strcmp(event.Key,'i') % visited port 7
%    assignin("caller","activeport",8)
% elseif strcmp(event.Key,'a') % visited port 8
%     assignin("caller","activeport",9)
% elseif strcmp(event.Key,'s') % visited port 9
%     assignin("caller","activeport",10)
% elseif strcmp(event.Key,'d') % visited port 10
%    assignin("caller","activeport",11)
% elseif strcmp(event.Key,'f') % visited port 11
%     assignin("caller","activeport",12)
% elseif strcmp(event.Key,'g') % visited port 12
%   assignin("caller","activeport",13)
% elseif strcmp(event.Key,'h') % visited port 13
%     assignin("caller","activeport",14)
% elseif strcmp(event.Key,'j') % visited port 14
%     assignin("caller","activeport",15)
% elseif strcmp(event.Key,'k') % visited port 15
%     assignin("caller","activeport",16)
elseif strcmp(event.Key,'m') % to manually trigger bell. automatically does on port switch
    bellSoundPath = 'sounds/service-bell_daniel_simion.mp3';
    playLocalAudioFile(bellSoundPath)  
elseif strcmp(event.Key,'return') % to define cleanup func, do this THEN ctrl c !!! IMPORTANT!!
     assignin("caller","closeupshop",1)
     print ('initiating cleanup func')
end

end
