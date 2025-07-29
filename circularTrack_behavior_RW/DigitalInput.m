classdef DigitalInput
    %DIGITALINPUT Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        Device
        NumberOfChannels
    end
    
    methods
        function obj = DigitalInput()
            d = daqlist("ni");
            deviceInfo = d{1, "DeviceInfo"};
            dq = daq("ni");
            obj.Device=dq;           
            DigitalIO=deviceInfo.Subsystems(3);
            chans=[6 2 7 3 8 4 20 19 23 21 18 17 16 15 24 14]; % input sensor ch 0 - 15 reference  on pcb
            for ichan=chans
                cn=DigitalIO.ChannelNames{ichan};
                addinput(dq, deviceInfo.ID, cn, "Digital");
            end

            %assigning values to class properties 
            obj.NumberOfChannels=numel(chans);
        end
        
        function d = readNI(obj,num)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            d=read(obj.Device,1,"OutputFormat","Matrix");
            if exist("num","var")
                d=d(num);
            end
        end
    end
end

