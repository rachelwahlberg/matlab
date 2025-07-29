% [PlaceMap, OccupancyMap] = PFClassic(Pos, SpkCnt, Smooth, nGrid, TopRate)
% Place field calculation "classic style" where a smoothed spike count map
% is divided by a smoothed occupancy map
%
% Pos is a nx2 array giving position in each epoch.  It should be in
% the range 0 to 1.  SpkCnt gives the number of spikes in each epoch.
%
% Smooth is the width of the Gaussian smoother to use (in 0...1 units).
%
% nGrid gives the grid spacing to evaluate on (should be larger than
% 1/Smooth). Because it's 2D should be three vals, [val1 val2 val3] for each
% dimension and then time.
%
% TopRate is for the maximum firing rate on the color map (if you display it)
% if you don't specify this it will be the maximum value
%
% optional output OccupancyMap is a smoothed occupancy map

function [pc3,TimeSpent,nSpikes,sTimeSpent,snSpikes,pc1,pc2] = PFClassic(Pos, SpkCount, Smooth, nGrid, Tbin, TopRate)
%function [pc3, nSpikes, sTimeSpent] = PFClassic(Pos, SpkCount, Smooth, nGrid, Tbin, TopRate)
% when the rat spends less than this many time points near a place,
% start to dim the place field
TimeThresh = .2;

% integrized Pos (in the range 1...nGrid
iPos = 1+floor(nGrid.*Pos/(1+eps));

% make unsmoothed arrays
TimeSpent = full(sparse(iPos(:,1), iPos(:,2),1, nGrid(1), nGrid(2)))*Tbin;
nSpikes = full(sparse(iPos(:,1), iPos(:,2), SpkCount, nGrid(1), nGrid(2)));

% do the smoothing
nGridmin=min(nGrid);
r = (-nGridmin:nGridmin)/nGridmin;
Smoother = exp(-r.^2/Smooth^2/2);
Smoother = Smoother/sum(Smoother);

sTimeSpent = conv2(Smoother, Smoother, TimeSpent, 'same');
snSpikes = conv2(Smoother, Smoother, nSpikes, 'same');


%SMOOTHING VERSION 1
pc1 = snSpikes./(sTimeSpent + eps); % version 1

%SMOOTHING VERSION 2
rr = 5;
ff = neuro.placeField.createSmoothingFilter(rr);
smoothedOcc = filter2(ff,TimeSpent,'same');
smoothedspike = filter2(ff,nSpikes,'same');
pc2 = smoothedspike./(smoothedOcc+eps);
%pc2(pc2 ==0) = nan;
%smoothedspike(obj.Empties) = nan;

%SMOOTHING VERSION 3
% smooth AFTER dividing raw spikes by smoothed occupancy
pc3 = conv2(Smoother, Smoother, nSpikes./(sTimeSpent + eps), 'same');

 

% NB not regularized to be mean firing rate in non-visited areas.

%if nargout==0
   % FireRate = PlaceMap*20000/512;
    if nargin<6
        TopRate = [];
    end
  %  TopRate = max(max(pc1)); %max out the saturation
    Gamma = 1; %higher gamma = more focus on saturated parts
    external.Placefields.PFPlot(pc1, sTimeSpent, TopRate,TimeThresh,Gamma);

    %end
%     if 0
%         colormap(gca, jet);
%         imagesc(PlaceMap);
%     else
%         FireRate = PlaceMap*20000/512;
%         %CoV = snSpikes.^-.5;
%         if nargin<5 | isempty(TopRate)
%             TopRate = max(FireRate(find(sTimeSpent>TimeThresh)));
%             if TopRate<1, TopRate=1; end;
%         end
%         if isempty(TopRate) TopRate = max(FireRate(:)); end;
%     	Hsv(:,:,1) = (2/3) - (2/3)*clip(FireRate'/TopRate,0,1);
%     	%Hsv(:,:,1) = (2/3) - (2/3)*FireRate'/MaxFireRate;
%     	%Hsv(:,:,3) = (sTimeSpent'/(max(sTimeSpent(:))+eps)).^.35;
%         %Hsv(:,:,3) = 1./(1+CoV');
%         Hsv(:,:,3) = 1./(1+TimeThresh./sTimeSpent');
%     	Hsv(:,:,2) = ones(size(FireRate'));
%     	image(hsv2rgb(Hsv));
%
%         % most annoying bit is colorbar
%         h = gca;
%         h2 = SideBar;
%         BarHsv(:,:,1) = (2/3) - (2/3)*(0:.01:1)';
%         BarHsv(:,:,2) = ones(101,1);
%         BarHsv(:,:,3) = ones(101,1);
%         image(0,(0:.01:1)*TopRate, hsv2rgb(BarHsv));
%         set(gca, 'ydir', 'normal');
%         set(gca, 'xtick', []);
%         set(gca, 'yaxislocation', 'right');
%         axes(h);
% %        keyboard
%     end
% end
