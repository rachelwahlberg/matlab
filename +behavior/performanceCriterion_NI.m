%%  Target pairs + general params

rat = 'Harry';

% sensor names from the PCB BOARD/PHYSICAL SENSORS, [# #] 
sensorPair = [2 3]; % 
%sensorPair = []; %PREFILL IN FOR EASY CHANGE MID SESSION
valvePair = [];
%solenoidPair = []; %PREFILL IN FOR EASY CHANGE MID SESSION

% Convert sensor name. send in as double, comes out as cell
sensors = convertSensorName(sensorPair);
valves = convertSensorName(valvePair);

%% receiving input from national instruments
% To read inputs from 16 sensors on NI board and
% output to up to 4 solenoids

di = DigitalInput();
do = DigitalOutput();

%% Switching rules

%initialize
sensor1 = sensors(1);
sensor2 = sensors(2); 
valve1 = valves(1);
valve2 = valves(2);

thresh = .85; % 85% performance

starttime = [];
if isempty(starttime); starttime = datetime; end

sensorTrig = 0;
flag = 0;

pumpOpen = 120; % in milliseconds

%%%%%%%%%%% If sensor1 is activated and it hasn't been triggered previously
if %%% sensor1 is activated
    %%% open valve 1 for pumpOpen length
    %%% record this happened
%     digitalWrite(valve2, HIGH);   // turn the LED on (HIGH is the voltage level)
%     digitalWrite(OpenEphys, HIGH);
%     delay(pumpOpen);               // wait for a second
%     digitalWrite(valve2, LOW);    // turn the LED off by making the voltage LOW 
%     digitalWrite(OpenEphys, LOW);
%     flag = 2;
%     SensorConf = 1;
%     Serial.println(SensorConf);
%     delay(pumpOpen);
%     sensorTrig = 1;
end

%%%%%%%%%%% If sensor1 is activated and it hasn't been triggered previously
if %%% sensor1 is activated
    %%% open valve 1 for pumpOpen length
    %%% record this happened
%     digitalWrite(valve2, HIGH);   // turn the LED on (HIGH is the voltage level)
%     digitalWrite(OpenEphys, HIGH);
%     delay(pumpOpen);               // wait for a second
%     digitalWrite(valve2, LOW);    // turn the LED off by making the voltage LOW 
%     digitalWrite(OpenEphys, LOW);
%     flag = 2;
%     SensorConf = 1;
%     Serial.println(SensorConf);
%     delay(pumpOpen);
%     sensorTrig = 1;
end


if %%% sensor2 is activated
    %%% open valve 2 for pumpOpen length
    %%% record this happened
end



continoo = 1;
n = 0;
totalcorrect = 0;











%% performance criterion




figure
hold on
title(['Performance de ' rat])
xlabel('minutes from start')
ylabel('performance')
while continoo == 1
    prompt = "Correct? 1 for yes, 0 for no";
    outcome = logical(input(prompt)); 
    totalcorrect = totalcorrect + outcome; 
    n=n+1;
    performance = totalcorrect/n;
    d = datetime;
    scatter(minutes(d-starttime),performance)
    disp(['Performance = ' num2str(performance)])

    
    if performance >= 0.85 
       disp('YOU HAVE HIT THRESHOLD !! CELEBRATION !!')
       %saveas(gcf,[char(datetime("today")) '_' rat '_performance.fig'])
    end
end

% 1 = white, outside, correct


    