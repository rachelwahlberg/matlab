classdef TimeIntervalZT < time.TimeInterval
    %TIMEINTERVALZT Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        ZeitgeberTime
    end
    
    methods
        function obj = TimeIntervalZT(varargin)
            %TIMEINTERVALZT Construct an instance of this class
            %   Detailed explanation goes here
            if isa(varargin{1},'time.TimeInterval')
                ti=varargin{1};
                startTime=ti.StartTimeAbs;
                sampleRate=ti.SampleRate;
                numberOfPoints=ti.NumberOfPoints;
                zt=varargin{2};
            elseif  isa(varargin{1},'time.TimeIntervalCombined')
                ti=varargin{1};
                startTime=ti.getStartTimeAbs;
                sampleRate=ti.getSampleRate;
                numberOfPoints=ti.getNumberOfPoints;
                zt=varargin{1}.getZ
            else
                startTime=varargin{1};
                sampleRate=varargin{2};
                numberOfPoints=varargin{3};
                zt=varargin{4};
            end
            obj@time.TimeInterval(startTime, sampleRate, numberOfPoints)
            if isduration(zt)
                obj.ZeitgeberTime= zt+obj.getDate;
            elseif isdatetime(zt)
                obj.ZeitgeberTime= zt;
            end
        end
        function S=getStruct(obj)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            S=getStruct@time.TimeInterval(obj);
            S.ZeitgeberTime=obj.ZeitgeberTime;
        end        
        function str=tostring(obj)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            str1=tostring@time.TimeInterval(obj);
            str2=char(obj.ZeitgeberTime,'yyyy-MM-dd--haa');
            str=sprintf('\tZT:%s %s',str2,str1);
        end

        function zt=getZeitgeberTime(obj)
            zt=obj.ZeitgeberTime;
        end
        function tps=getTimePointsZT(obj)
            st=obj.getStartTimeZT;
            en=obj.getEndTimeZT;
            tps=linspace(st,en,obj.NumberOfPoints);
        end
        function st=getStartTimeZT(obj)%gives a duration object
            st1=obj.getStartTimeAbs;
            zt=obj.getZeitgeberTime;
            st=st1-zt;
        end
        function zttimes=getZTTimeForSamples(obj,samples)
            %same function as "getRealTimeFor" but subtracts the zt time
            %(for example, noon, when the dark/light transition is) from
            %the real time stamps to output zt times.
            %             if numel(samples)>1e5
            %                 tps=obj.getTimePointsInAbsoluteTimes;
            %                 zttimes=tps(samples);
            %             else
            zt = obj.getZeitgeberTime;
            for isample=1:numel(samples)
                sample=samples(isample);
                zttimes(isample) = obj.getRealTimeFor(sample)-zt;
            end
            %             end
        end
    end
end

