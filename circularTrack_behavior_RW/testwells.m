function [] = testwells(pcb_valve,pumpOpen)

% to trigger wells without needing to activate with a port
%pcb_valve = 3; % 0, 1, 2, or 3
%pumpOpen =  0.5; % in seconds, length of time water is delivered (small because of the gravity system)
%don't go lower than 0.02! solenoid won't work. sometimes unreliable up to
%0.7 as well.

valve = convertSensorName(pcb_valve);

do = DigitalOutput(); %out to up to 4 solenoids
do.toggleNI(valve,pumpOpen) % output,time(seconds)



