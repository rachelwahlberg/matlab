%pc_velocityfiltercomparison
%to run from methodsFigure (in Pipelines)

%plots spikelocations, occupancy map, spike map, place map for all
%velocities (top row) and velocities > 10cm (bottom row)


%%%%% uses PFClassic: %%%%%%%
%pos is nx2 array, NORMALIZED between 0 and 100.
%spksBinned gives the number of spikes in each epoch.
%Smooth is width of gaussian smoother (in 0 ... 1 units)
% nGrid gives grid spacing (should be larger than 1/smooth)

import position.* %otherwise when I call position.getDirectionsIdx it's calling a different position func?

%%%%%% default vals %%%%%%%%
timeWindow = [obj.positionData.time.getStartTimeAbs obj.positionData.time.getEndTimeAbs];
smoothVal = .02; %out of 100
nGrid = [100 100]; %one val for each dimension + for time
timebinsize = .01;% tbins in units of seconds, e.g. 10 ms.

%%%%%%%%% time %%%%%%%%%%%%%%%%
positionTimeZT = obj.positionData.time.getTimePointsZT;

%%%%%%%% Normalize position  %%%%%%%%%%%%%%
x = obj.positionData.data.X;
z = obj.positionData.data.Z;
x = x - min(x); %make all values positive, aligned to zero
z = z - min(z);
%
normX = x/max(x); %normalize to make all vals between 0 and 1.
normZ = z/max(z);

%%%%%%%%%% Get binned position %%%%%%%%%%%
%t = obj.positionData.time.getTimePointsInSamples/posSR; %in seconds
tbegin = positionTimeZT(1); tend = positionTimeZT(end);
nBins = round(seconds(tend-tbegin)/timebinsize);
xzt(:,1) = normX;
xzt(:,2) = normZ;
xzt(:,3) = seconds(positionTimeZT);

positionBinned = external.Placefields.binpos(xzt,nBins);
histedges = positionBinned(:,3)-timebinsize;
histedges(end+1) = histedges(end)+timebinsize;

%%%%%%%%% filter by... PART 2 %%%%%%%%%%%

%%% theta
%  belowthreshtheta = filterbytheta();

%%%%%% spiking info %%%%%%%%%%
sf = neuro.spike.SpikeFactory.instance(); %Spike Factory (Kaya code)
sa = sf.getPhyOutputFolder(obj.phypath); %Spike Array (Kaya code)
su = sa.getSpikeUnits; % pull spiketimes per unit.
nUnits = length(su);

try
    classTypes = unique(sa.ClusterInfo.class);
catch
    classTypes = 1;
end

nclassTypes = length(classTypes);

%%%%%% for each unit, plot place field %%%%%%%
figure
grid_height = 2; grid_width = 4;
for c = 1:nclassTypes
    class = classTypes(c);
    for unit = 1:nUnits

        %%%%%% get spikeCount %%%%%%%%
        % spsBinned is binned with same bins as binnedPos above

        spksZT = seconds(su(unit).getTimesZT)';  %spikes in ZT time.

        spksBinned = histcounts(spksZT,histedges)';
        spksBinned(end)=0; spksBinned(1)=0; % find number of spikes per unit time.

        spkslogical = logical(spksBinned);
        xprevelocityfilter = positionBinned(spkslogical,1);
        zprevelocityfilter = positionBinned(spkslogical,2);

        spksBinned_prefilter = spksBinned;

        %%%%%%%%% % velocityfilter = external.velocityfilter(binnedPos)
        velocityThresh = 10; %in 10 cm/hour
        [~,~,spksBinned] = position.velocityfilter_2D(positionBinned,velocityThresh,spksBinned); %in cm/s

        spkslogical = logical(spksBinned);
        xpostvelocityfilter = positionBinned(spkslogical,1);
        zpostvelocityfilter = positionBinned(spkslogical,2);

        h1 = tiledlayout(grid_height,grid_width);
        title(h1,['Unit ' num2str(unit)]);
        pos = 1; h = 1; w = 1; ax1 = nexttile(pos,[h,w]);
        pos = 2; h = 1; w = 1; ax2 = nexttile(pos,[h,w]);
        pos = 3; h = 1; w = 1; ax3 = nexttile(pos,[h,w]);
        pos = 4; h = 1; w = 1; ax4 = nexttile(pos,[h,w]);
        pos = 5; h = 1; w = 1; ax5 = nexttile(pos,[h,w]);
        pos = 6; h = 1; w = 1; ax6 = nexttile(pos,[h,w]);
        pos = 7; h = 1; w = 1; ax7 = nexttile(pos,[h,w]);
        pos = 8; h = 1; w = 1; ax8 = nexttile(pos,[h,w]);

        % AX1: plots position data vertically over time, and then the locations of unit firing
        axes(ax1)
        plottitle = 'Spk locs. (all velocities)';
        plot.plotPlaceFields('spkLocations',plottitle,x,z,xprevelocityfilter,zprevelocityfilter)

        % AX3: Plots occupancy and firing rate map
        axes(ax4)
        hold on

        %Get and plot place field %%%%%%%
        [~,~,nSpikes,sTimeSpent,snSpikes] = ...
            external.Placefields.PFClassic(positionBinned(:,1:2),spksBinned_prefilter,smoothVal,nGrid,timebinsize);
        title(['Place map, nspikes = ' num2str(sum(sum(nSpikes)))])
        % AX4:
        axes(ax2)
        plottitle = 'Smoothed Occ., all vel.';
        plot.plotPlaceFields('imagescMap',plottitle,sTimeSpent)

        % AX5:
        axes(ax3)
        plottitle = 'Smoothed Spikes,all vel.';
        plot.plotPlaceFields('imagescMap',plottitle,snSpikes)

        % AX2: plots position data vertically over time, and then the locations of unit firing
        axes(ax5)
        plottitle = 'Spk locs. (> 10 cm velocities)';
        plot.plotPlaceFields('spkLocations',plottitle,x,z,xpostvelocityfilter,zpostvelocityfilter)

        % AX3: Plots occupancy and firing rate map
        axes(ax8)
        hold on

        %Get and plot place field %%%%%%%
        [~,~,nSpikes,sTimeSpent,snSpikes] = ...
            external.Placefields.PFClassic(positionBinned(:,1:2),spksBinned,smoothVal,nGrid,timebinsize);
        title(['Place map, nspikes = ' num2str(sum(sum(nSpikes)))])
        % AX4:
        axes(ax6)
        plottitle = 'Smoothed Occ., >10cm';
        plot.plotPlaceFields('imagescMap',plottitle,sTimeSpent)

        % AX5:
        axes(ax7)
        plottitle = 'Smoothed Spikes, >10cm';
        plot.plotPlaceFields('imagescMap',plottitle,snSpikes)

        pause
        clf
    end
end

