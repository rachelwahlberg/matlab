function [clockwiseIdx,counterIdx] = getDirectionsIdx(Turns,spksBinned)
% send in Turns, which is created in directionsCircularTrackRadians. 
% currently being called from methodsFigure. 231221

clockwiseIdx = 0;
counterIdx = 0;
for t = 1:height(Turns)
    if t == height(Turns)
        break
    end

    if Turns.Direction(t) == 1 %clockwise
    tmpIdx = Turns.Index(t):(Turns.Index(t+1)-1);
    clockwiseIdx = [clockwiseIdx;tmpIdx'];
    end

    if Turns.Direction(t) == 0 %counterclockwise
    tmpIdx = Turns.Index(t):(Turns.Index(t+1)-1);
    counterIdx = [counterIdx;tmpIdx'];
    end
end

clockwiseIdx(1) = [];%get rid of the zeros
counterIdx(1) = [];

if exist("spksBinned","var")







end
        





























end

