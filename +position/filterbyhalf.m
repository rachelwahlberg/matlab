function [newZT,newX,newZ] = filterbyhalf(positionTimeZT,x,z,half)

if strcmp(half,"Both") == 1
    newZT = positionTimeZT;
    newX = x;
    newZ = z;
elseif strcmp(half,"First") == 1

    firsthalfindex = 1:floor(length(positionTimeZT)/2);
    newZT = positionTimeZT(firsthalfindex);

    newX = x(firsthalfindex);
    newZ = z(firsthalfindex);

elseif strcmp(half,"Second") == 1


    firsthalfindex = 1:floor(length(positionTimeZT)/2);
    secondhalfindex = (length(firsthalfindex)+1):length(positionTimeZT);
    newZT = positionTimeZT(secondhalfindex);

    newX = x(secondhalfindex);
    newZ = z(secondhalfindex);
end

end


