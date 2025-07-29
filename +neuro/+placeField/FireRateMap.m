classdef FireRateMap < neuro.placeField.OccupancyMap
    %FIRERATEMAP Summary of this class goes here
    %   Detailed explanation goes here

    properties
       % SpikePositions
        OccupancyMap
        OccupancyMapSmooth
        SpikeUnitTracked
        SmoothedSpike
    end

    methods
        %function obj = FireRateMap(occupancyMap,spikePositions) %constructor
            function obj = FireRateMap(occupancyMap,spiketimesZT) %constructor
            %lay spiking activity over occupancy map

            if nargin>0
                fnames=fieldnames(occupancyMap);
                for ifn=1:numel(fnames)
                    obj.(fnames{ifn})=occupancyMap.(fnames{ifn});
                end
                %obj.SpikePositions=spikePositions;
                obj.PositionData=occupancyMap.PositionData;
                obj.OccupancyMap=occupancyMap.MapOriginal;
                obj.OccupancyMapSmooth=occupancyMap.MapSmooth;
                if ~isempty(obj.XEdges)|| ~isempty(obj.ZEdges)
                    obj=obj.getSpikeOccupancy(spiketimesZT,obj.XEdges,obj.ZEdges);
                else
                    obj=obj.getSpikeOccupancy;
                end
            end
        end

        function [] = plot(obj)
            h = pcolor(obj.MapOriginal');
       %     axis equal
      %      colorbar
            set(h,'linestyle','none')
       %     axis equal
        end

        function [] = plotSmooth(obj)
            h = pcolor(obj.MapSmooth');
       %     axis equal %%% BRING THESE BACK
       %     colorbar
            set(h,'linestyle','none')
 %           axis equal
        end

        function [] = plotUK(obj) %kaya version
            ms=obj.OccupancyMap;
            ms(ms<eps)=0;
            alpha1=log(ms); 
            x=[min(obj.PositionData.data.X) max(obj.PositionData.data.X)];
            z=[min(obj.PositionData.data.Z) max(obj.PositionData.data.Z)];
            imagesc(x,z,obj.MapOriginal,AlphaDataMapping="scaled",AlphaData=alpha1);  
        end
        
        function [] = plotSmoothUK(obj)
            ms=obj.OccupancyMap;
            ms(ms<eps)=0;
            alpha1=log(ms);
            x=[min(obj.PositionData.data.X) max(obj.PositionData.data.X)];
            y=[min(obj.PositionData.data.Z) max(obj.PositionData.data.Z)];
            imagesc(x,y,obj.MapSmooth,AlphaDataMapping="scaled",AlphaData=alpha1);
        end

        function obj = calculate(obj,xedges,zedges)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            obj= calculate@neuro.placeField.OccupancyMap(obj);
            spikePositions=obj.SpikePositions;
            pd=obj.PositionData.data;
            PosSpk=[spikePositions.X spikePositions.Z];
            if nargin==1
                nGrid(1)=round((max(pd.X)-min(pd.X))/obj.SpatialBinSizeCm);
                nGrid(2)=round((max(pd.Z)-min(pd.Z))/obj.SpatialBinSizeCm);
                [obj.MapOriginal]=histcounts2( ...
                    PosSpk(:,1),PosSpk(:,2),nGrid);
            else
                [obj.MapOriginal]=histcounts2( ...
                    PosSpk(:,1),PosSpk(:,2),xedges,zedges);
            end
            obj.MapOriginal=obj.MapOriginal'./seconds(obj.TimeBin);

            % do the smoothing
            obj.MapSmooth=imgaussfilt(obj.MapOriginal,obj.Smooth);
        end

        function obj = getSpikeOccupancy(obj,spiketimesZT,xedges,zedges)

            %%%%%%%%% pre allocate %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            alphas = nan(size(obj.OccupancyMap));%calculated in OccupancyMap.calculateAJ
            betas = nan(size(obj.OccupancyMap));
            FRmaps = nan(size(obj.OccupancyMap));


            %%%%%%% get bins (same as in occupancy map creation) %%%%%%%%%%
            pd=obj.PositionData.data;
            %sr = obj.PositionData.time.getSampleRate;
            x = pd.X; z = pd.Z;
            %Pos1=[pd.X pd.Z];

            binsize = obj.SpatialBinSizeCm;
            if nargin==1 % no xedges/zedges provided
                xedges = min(x):binsize:max(x);      % The x-coordinate bins for the grid
                zedges = min(z):binsize:max(z);      % The y-coordinate bins for the grid
                xedges = xedges(1:end-1) + 0.5*binsize;
                zedges = zedges(1:end-1) + 0.5*binsize;
            end
            nxedges = length(xedges);
            nzedges = length(zedges);
            xedgesInd = 1:nxedges;
            zedgesInd = 1:nzedges;

            xb = interp1(xedges,xedgesInd,x,'nearest');
            zb = interp1(zedges,zedgesInd,z,'nearest');

            %%%%%%%% Get spike occupancy %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            rr = 3;
            ff = neuro.placeField.createSmoothingFilter(rr);

            behaviortime = obj.PositionData.time.getTimePointsZT;
            xint = interp1(behaviortime,xb,spiketimesZT,'nearest');   % Align the position of the animal when the cell spikes to an x-grid position
            zint = interp1(behaviortime,zb,spiketimesZT,'nearest');   % Align the position of the animal when the cell spikes to an y-grid position

            spikeocc = zeros(nxedges,nzedges);
            for iX = 1:nxedges
                % Find all X positions in the bin
                tmpX = (xint == iX);
                for iZ = 1:nzedges
                    % Find all Y positions in the bin
                    tmpZ = (zint == iZ);
                    spikeocc(iX,iZ) = sum(tmpX.*tmpZ);  % Since tmpX and tmpY are binary, multiplying them will yield 1 only if both are 1.  This is effectively an AND function.
                end
            end
            smoothedspike = filter2(ff,spikeocc,'same');
          %  smoothedspike(obj.Empties) = nan;

            obj.MapOriginal=spikeocc./obj.OccupancyMap;
            obj.MapSmooth=smoothedspike./obj.OccupancyMapSmooth;% This is your standard spikes per unit time based map.  It's not statistically correct.  But it's close and easy.
            obj.MapOriginal(obj.Empties) = nan;
            obj.MapSmooth(obj.Empties) = nan;
            obj.SmoothedSpike = smoothedspike;
        end

        function pfm = getPlaceFieldMap(obj)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            pfm1=neuro.placeField.PlaceFieldMap(obj);
            pfm=pfm1.getPlaceFieldMapMeasures;
        end
    end
end

