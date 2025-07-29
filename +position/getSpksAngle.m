function [spksAngle,spksT] = getSpksAngle(spksBinned,positionBinned)

% feed in spksBinned for which bins to get, and then positionBinned for xy coordinates. 

index = logical(spksBinned);
spksX = positionBinned(index,1);
spksY = positionBinned(index,2);
spksT = positionBinned(index,3);

spksX = spksX - max(spksX)/2;
spksY = spksY - max(spksY)/2;

spksAngle = atan2(spksY,spksX);






