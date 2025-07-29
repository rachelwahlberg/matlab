function [] = filterbytheta



    % get eeg states
            % eegStates = getStates(obj,timeWindow);
            basepath = [obj.basepath '/'];
            % rejectCh = [97:127];
            % thetaCh = [64:74];
            % SleepState = SleepScoreMaster_RW(obj.basepath,'Notch60Hz',1,'overwrite',true);%,'rejectChannels',rejectCh,'thetaChannels',thetaCh,'overwrite',
            % Get theta power at each point in

            [thetaPow,thetaTfm] = getPowerTimeWindow(obj,timeWindow,[150 200]);
            lfpTimeInterval = obj.LFP.TimeIntervalCombined.getTimeIntervalForTimes(timeWindow);
            lfpTime = seconds(lfpTimeInterval.getTimePointsZT);

            thetaBinned = interp1(lfpTime,thetaPow.Values,positionBinned(:,3)); %find theta power at the times where lfpTime = binnedPos time.

            thetaThresh = 1;
            zTheta = (thetaBinned- nanmean(thetaBinned))/nanstd(thetaBinned);
            %figure
            %thetaHist = histogram(zTheta,100);
            belowthreshTheta = find(zTheta < 0); %ARBITRARY IN MIDDLE RIGHT NOW, no bimodality