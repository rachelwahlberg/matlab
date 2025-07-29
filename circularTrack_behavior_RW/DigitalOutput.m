classdef DigitalOutput
    %IOOUTPUT Summary of this class goes here
    %   Detailed explanation goes here
    
    %define properties of class
    properties
        VoltageRange
        SamplingRate
        NumberOfChannels
        Device
    end
    
    %deine methods of class - functions that implement class methods
    methods

        %sets the scanning rate of a class instance
        function obj = DigitalOutput(Rate)
            d = daqlist("ni");
            deviceInfo = d{1, "DeviceInfo"};
            dq = daq("ni");
            if ~exist('Rate','var')
                dq.Rate = 500000;
            else
                dq.Rate=Rate;
            end
            obj.Device=dq;           
            DigitalIO=deviceInfo.Subsystems(3);
            chans=[10 9 12 11]; % in 0 1 2 3 order
            for ichan=chans
                cn=DigitalIO.ChannelNames{ichan};
                addoutput(dq, deviceInfo.ID, cn, "Digital");
            end

            %assigning values to class properties 
            obj.NumberOfChannels=numel(chans);
            obj.SamplingRate=dq.Rate;
        end
        
        %writes 
        function [] = write(obj,points)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            write(obj.Device,points);
        end
        function [] = toggleNI(obj,eventPin,dur)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            if ~exist('dur','var')
                dur=.100;%20 ms
            end
            valve=[0 0 0 0];
            obj.Device.write(valve);
            valve(eventPin)=1;
            obj.Device.write(valve);
            pause(dur);
            valve(eventPin)=0;
            obj.Device.write(valve);
        end
    end
end

