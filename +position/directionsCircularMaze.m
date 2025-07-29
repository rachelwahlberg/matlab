function [clockwiseTimes,counterTimes] = directionsCircularMaze(positionData,rawData)

%Send in PositionData object. 
% if second is sent in, it's for the speed calculation (for ex, if
% linearized positionData is var1 and you want non linearized for speed
% calculation, you send that in as a second argument)

%   This function takes position data and partitions the maze into
%   sections.
%
%   INPUTS:
%      positionData object.
%
x = positionData.data.X;
z = positionData.data.Z;
smoothingwindowinseconds = 1;
speed = rawData.getSpeed(smoothingwindowinseconds);

x(speed.Values<10) = nan;
z(speed.Values<10) = nan;

%% Get xy coordinate bounds for maze sections.
zmax = max(z); xmin = min(x);
xmax = max(x); zmin = min(z);
%zmid = 0.5*(zmin + zmax);
xmid = 0.5*(xmin + xmax);

clockwise = logical(zeros(size(x)));
counter = logical(zeros(size(x)));
for p = 1:length(x)
    if isnan(x(p))
        continue
    end

    if p == length(x) %last point
        if istrue(counter(p-1))
            counter(p) = true;
        elseif istrue(clockwise(p-1))
            clockwise(p) = true;
        end
        break
    end


    if x(p)>=xmid & z(p+1)>z(p)
        counter(p) = true; %going CLOCKWISE
    elseif x(p)>=xmid & z(p+1)<z(p)
        clockwise(p) = true; %going COUNTERCLOCKWISE
    elseif x(p)<xmid & z(p+1)>z(p)
        clockwise(p) = true; %going COUNTERCLOCKWISE
    elseif x(p)<xmid & z(p+1)<z(p)
        counter(p) = true; %going CLOCKWISE
    end

end

zt = positionData.time.getTimePointsZT;
clockwiseTimes = zt(clockwise);
counterTimes=zt(counter);


%get the timestamps for the clockwise portion
