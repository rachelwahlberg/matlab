function [dirLogical,clockwise_spksBinned,counter_spksBinned] = counterClockwise_positionBinned(positionBinned,Turns,spikesBinned)
%must be AT LEAST the threshold - eliminates slower values

%take the differential of x values and then y values to determine points
%for which speed is below the provided threshold
binnedTime = positionBinned(:,3);
turnTimestamps = Turns.TimestampZT; %includes counter + clockwise
turnDir = Turns.Direction; % 1 is clockwise, 0 is counter

dirLogical = nan(length(positionBinned),1); %should have no nans at end
for t = 1:length(turnTimestamps)

    [~,idx] = min(abs(binnedTime-turnTimestamps(t))); %find closest timestamp in positionBinned

    if t ~= length(turnTimestamps)
        [~,nextIdx] = min(abs(binnedTime-turnTimestamps(t+1))); %find closest timestamp in positionBinned

        if turnDir(t) == 0 % if turn is counterclockwise
            dirLogical(idx:nextIdx-1) = 0; %next segment is counter
        else %if turnDir(t) == 1 % if turn is clockwise
            dirLogical(idx:nextIdx-1) = 1; %next segment is clockwise
        end

        if t == 1
            if turnDir(t) == 0 % if first turn is counterclockwise
                dirLogical(1:idx-1) = 1; %before then was clockwise
            else% if turnDir(t) == 1 % if first turn is clockwise
                dirLogical(1:idx-1) = 0; %before then was counterclockwise
            end
        end

    else %if t == length(turnTimestamps)
        if turnDir(t) == 0 % if last turn is counterclockwise
            dirLogical(idx:end) = 0; %rest is counterclockwise
        else %if turnDir(t) == 1 % if last turn is clockwise
            dirLogical(idx:end) = 1; %rest is clockwise
        end
    end
end

dirLogical = logical(dirLogical);

if exist("spikesBinned","var") %get out spike bins that are clockwise/counter
    clockwise_spksBinned = spikesBinned;
    counter_spksBinned = spikesBinned;

    clockwise_spksBinned(~dirLogical) = 0; %any counter bins become xero
    counter_spksBinned(dirLogical)=0; % any clockwise bins become zero
end

end