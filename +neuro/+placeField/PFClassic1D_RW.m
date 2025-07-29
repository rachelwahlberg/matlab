% [PlaceMap, OccupancyMap] = PFClassic(Pos, SpkCnt, Smooth, nGrid, TopRate)
% Place field calculation "classic style" where a smoothed spike count map
% is divided by a smoothed occupancy map
%
% Pos is a nx2 array giving position in each epoch.  It should be in
% the range 0 to 1.  SpkCnt gives the number of spikes in each epoch.
%
% Smooth is the width of the Gaussian smoother to use (in 0...1 units).
%
% nGrid gives the grid spacing to evaluate on (should be larger than 1/Smooth)
%
% TopRate is for the maximum firing rate on the color map (if you display it)
% if you don't specify this it will be the maximum value
%
% optional output OccupancyMap is a smoothed occupancy map

function [PlaceMap, spikeMap, TimeSpent] = PFClassic1D_RW(Pos, SpkCount, Smooth, nGrid,Tbin,plotplacefield)
% for circular track - vals come in -pi to pi.

    % when the rat spends less than this many time points near a place,
    % start to dim the place field
    TimeThresh = .20; %Rw added
    Pos(Pos<0) = Pos(Pos<0)+(2*pi);%to convert from -pi to pi range to 0 to 2pi 
    Pos = Pos/max(Pos); %not subtracting min first on purpose, because I want to keep true zero

    % integrized Pos (in the range 1...nGrid
    iPos = 1+floor(nGrid*Pos/(1+eps));

    % make unsmoothed arrays
    TimeSpent = full(sparse(ones(size(iPos)),iPos, 1, nGrid, nGrid));
    spikeMap = full(sparse(ones(size(iPos)),iPos, SpkCount, nGrid, nGrid));

    TimeSpent = TimeSpent(1,:)*Tbin;
    spikeMap = spikeMap(1,:);
    
    % do the smoothing
    r = (-nGrid:nGrid)/nGrid;
    Smoother = exp(-r.^2/Smooth^2/2);
    Smoother = Smoother/sum(Smoother);
    
    %sTimeSpent = conv(TimeSpent,Smoother,'same');
    %spikeMap = conv(nSpikes, Smoother,'same');

    PlaceMap = spikeMap./(TimeSpent+eps);
    PlaceMap = conv(PlaceMap, Smoother,'same');

    % NB not regularized to be mean firing rate in non-visited areas.

    if plotplacefield == 1
        if nargin<7
            TopRate = [];
        end
        Gamma = 1; %higher gamma = more focus on saturated parts
        external.Placefields.PFPlot(PlaceMap,TimeSpent,TopRate,TimeThresh,Gamma);
    end

end