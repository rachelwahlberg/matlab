classdef SpikeUnitTracked < neuro.spike.SpikeUnit
    %SPIKEUNITTRACKED Summary of this class goes here
    %   Detailed explanation goes here

    properties
        PositionData
        InfoPosition
    end

    methods
        function obj = SpikeUnitTracked(spikeunit,positionData)
            %SPIKEUNITTRACKED Construct an instance of this class
            %   Detailed explanation goes here
            if nargin>0
                fnames=fieldnames(spikeunit);
                for ifn=1:numel(fnames)
                    fn=fnames{ifn};
                    try
                        obj.(fn)=spikeunit.(fn);
                    catch ME
                        
                    end
                end
                obj.PositionData=positionData;
                obj.TimesInSamples=[spikeunit.TimesInSamples obj.getPositions];
            end
        end

        function [positions] = getPositions(obj)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes hereobj.
            pos=obj.PositionData;
            pt=seconds(pos.time.getTimePointsZT)';
            sts=seconds(obj.getTimesZT)'; %all spikes, not divided by label
           % sts2=seconds(obj.getTimesZT)';  

            pt(find(isnan(pos.data.X(:)))) = nan;
            
            % Create a KDTree object from pt
            tree = KDTreeSearcher(pt);
            % Find the indices of the closest values in pt for each value
            % in sts % could be way off because a lot of it isn't on the
            % track itself? 
            locinpos = knnsearch(tree, sts);
            threshold=max(1./[obj.Time.getSampleRate ...
                pos.time.getSampleRate]); %max between 1/behavior sr and 1/lfp sr
            locinpos(abs(sts - pt(locinpos)) > threshold) = NaN;
            positions=array2table(nan(length(sts),3), ...
                VariableNames=pos.data.Properties.VariableNames);
            % for isp=1:numel(sts)
            %     st=sts(isp);
            %     [g, idx]=min(abs(st-pt));
            %     if g<max(1./[obj.TimeIntervalCombined.getSampleRate ...
            %             pos.time.getSampleRate])
            %         locinpos(isp)=idx;
            %     end
            % end
            idxvalid=~isnan(locinpos);
            positions(idxvalid,:)=pos.data(locinpos(idxvalid),:);
        end

        function pd = getPositionData(obj)
            pd=obj.PositionData;
        end
        function [] = plotOnTimeTrack(obj,speedthreshold)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            alpha=.3;
            track=obj.PositionData;
            track.plot;hold on
            times=obj.getAbsoluteSpikeTimes;
            timesmin=minutes(times-track.time.getStartTimeAbs);
            timeratio=timesmin./...
                minutes(track.time.getEndTime-track.time.getStartTimeAbs);
            [X,Y,Z,idx]=track.getLocationForTimesBoth( ...
                times,speedthreshold);
            color=linspecer(11);
            try
                s1=scatter(Z(idx),timesmin(idx),50,color(round( ...
                    timeratio(idx)*10)+1,:),'filled');
                s2=scatter(Z(~idx),timesmin(~idx),5,[0 0 0],'filled');
            catch
                error
            end
            ax=gca;
            ax.YDir='reverse';
            legend off
            s1.MarkerFaceAlpha=alpha;
            s1.MarkerEdgeAlpha=alpha;
            s2.MarkerFaceAlpha=alpha/3;
            s2.MarkerEdgeAlpha=alpha/3;
            str=obj.addInfo(idx);
        end

        

        function [] = plotOnTrack2D(obj) %
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            numPointsInPlot=50000;
            positionData=obj.PositionData;
            positionData.plot2D(numPointsInPlot);hold on

            if numel(obj.TimesInSamples)>0
                [~,idx]=positionData.getPositionForTimes( ...
                    obj.getAbsoluteSpikeTimes);
                positionData.plot2DMark(idx);hold on
            end
        end
        function [] = plotSpikes2D(obj)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            positionData=obj.PositionData;
            if numel(obj.TimesInSamples)>0
                [~,idx]=positionData.getPositionForTimes( ...
                    obj.getAbsoluteSpikeTimes);
                positionData.plot2DMark(idx);hold on
            end
        end

        function [ax] = plotOnTrack3D(obj,events,color)
            %Plots evolution of position across time with time on the z
            %axis, 2D position on the xy axes. overlays events (such as
            %spiketimes)

            if ~exist("color","var")
                color=[];
            end

            numPointsInPlot=50000;
            track=obj.PositionData;
            track.plot3DtimeContinuous(numPointsInPlot);hold on

            if numel(obj.TimesInSamples)>0
                [~,idx]=track.getPositionForTimes( ...
                    obj.getAbsoluteSpikeTimes);
                track.plot3DMark(idx,color);hold on
            end
            ax=gca;
            ticd=obj.PositionData.time;
            t_org=seconds(ticd.getTimePoints);
            ax.ZLim=[0 abs(t_org(2)-t_org(end))];
            ax.ZDir="reverse";
            if exist("events","var") %  for a dot at the time of a specific event (similar idea to plotplane)
                [~,idx]=track.getPositionForTimes(events.ruleSwitch);
                track.plot3DMark(idx,'r')
            end
        end

        function [] = plotPlaneonTrack3D(obj,events,color)
            % plots a plane (with the purpose of being over
            % plotonTrack3D) at certain events (goal being to plot a
            % plane at the time of epoch change, for ex)
            positionData = obj.PositionData;
            if ~exist("color","var")
                color=[0.3010 0.7450 0.9330];
            end

            [~,idx]=positionData.getPositionForTimes(events); % can feed in as absolute time, zt time, etc
            positionData.plot3DEventPlane(idx,color);
        end
        
        function [frm] = getFireRateMap(obj,xedges,zedges)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            %track=obj.PositionData;
            if nargin>1
                om=obj.PositionData.getOccupancyMap(xedges,zedges);
            else
                om=obj.PositionData.getOccupancyMap();
            end
  
            %spiketimesZT = obj.getSpikeTimesZT;
            spiketimesZT = obj.getTimesZT;
            frm = neuro.placeField.FireRateMap(om,spiketimesZT);
            frm.SpikeUnitTracked=obj;

%             [sTimes,~]=track.getPositionForTimes( ...
%                 obj.getAbsoluteSpikeTimes);
%             
% frm=om+sTimes; 
% 


        end

        function [tall] = getSpikePositionTable(obj)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            track=obj.PositionData;
            
            [sTimes,~]=track.getPositionForTimes( ...
                obj.getAbsoluteSpikeTimes);
            t1=seconds(obj.getTimesZT)';
            t2=array2table(t1,VariableNames={'TimeZT'});
            tall=[sTimes t2];
        end


        function [cache,obj] = getPositionHeld(obj,cache)
            [cache,key]=cache.hold( ...
                obj.PositionData);
            obj.PositionData=key;

        end

        function st=getSpikeTimesZT(obj)
            ts=obj.getSpikeArrayWithAdjustedTimestamps; %currently does nothing?
            st1=seconds(double(obj.TimesInSamples.SpikeTimes)/ ...
                ts.Time.getSampleRate);
            st=obj.Time.getStartTimeZT+st1;
        end

         function obj = getSpikeArrayWithAdjustedTimestamps(obj)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            st=obj.TimesInSamples.SpikeTimes;
            %spiketimessample=st.SpikeTimes;
            ticd=obj.Time;
            adjustedspiketimessample= ...
                ticd.adjustTimestampsAsIfNotInterrupted(st);
            obj.TimesInSamples.SpikeTimes=adjustedspiketimessample;
            if isa(ticd,'time.TimeIntervalCombined')
                ti=ticd.timeIntervalList.get(1);
            elseif isa(ticd,'time.TimeInterval')
                ti=ticd;
            end
            ti.NumberOfPoints=ticd.getNumberOfPoints;
            obj.Time = ti;
        end

        function [obj] = getPositionReload(obj,cache)
            try
                obj.PositionData=cache.get(obj.PositionData);
            catch ME

            end

        end

    end
end

