classdef OccupancyMap
    %OCCUPANCYMAP Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        TimeBin
        SpatialBinSizeCm
        Smooth
        PositionData
        MapOriginal
        MapSmooth
        Units
        XEdges
        ZEdges
        Empties
    end
    
    methods
        function obj = OccupancyMap(positionData,srate,xedges,zedges)
            %OCCUPANCYMAP Construct an instance of this class
            %   Detailed explanation goes here

         
            if nargin>0
                
                obj.PositionData = positionData;
                obj.TimeBin=seconds(1/srate);
                obj.Smooth= 4; %5
                obj.SpatialBinSizeCm=2;
                if nargin>2
                    obj=obj.calculateAJ(xedges,zedges);
                else
                    obj=obj.calculateAJ();
                end
            end
        end

        function [] = plot(obj)
            ms=obj.MapSmooth;
            ms(ms<eps)=0;
            alpha1=log(ms);
            x=obj.getXLim;
            y=obj.getZLim;
            imagesc(x,y,obj.MapOriginal,AlphaDataMapping="scaled",AlphaData=alpha1);    
        end
        function [] = plotSmooth(obj)
            ms=obj.MapSmooth;
            ms(ms<eps)=0;
            alpha1=log(ms);
            x=obj.getXLim;
            y=obj.getZLim;
           imagesc(x,y,obj.MapSmooth,AlphaDataMapping="scaled",AlphaData=alpha1);    
        end
        function obj = setTimeBin(obj,val)
            if ~isduration(val)
                error('value should be in Duration')
            end
            obj.TimeBin=val;
            obj=obj.calculate;
        end
        function obj = setSmooth(obj,val)
            obj.Smooth=val;
            obj=obj.calculate;
        end
        function obj = setSpatialBinSizeCm(obj,val)
            obj.SpatialBinSizeCm=val;
            obj=obj.calculate;
        end
        
        function obj = calculate(obj,xedges,zedges)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            pd=obj.PositionData.data;
            Pos1=[pd.X pd.Z];
            if nargin==1
                nGrid(1)=round((max(pd.X)-min(pd.X))/obj.SpatialBinSizeCm);
                nGrid(2)=round((max(pd.Z)-min(pd.Z))/obj.SpatialBinSizeCm);
                try
                    [obj.MapOriginal,Xedges,Zedges]=histcounts2( ...
                    Pos1(:,1),Pos1(:,2),nGrid);
                catch ME
                    
                end
            else
                [obj.MapOriginal,Xedges,Zedges]=histcounts2( ...
                    Pos1(:,1),Pos1(:,2),xedges,zedges);
            end
            obj.MapOriginal=obj.MapOriginal';
            % do the smoothing
            obj.MapSmooth=imgaussfilt(obj.MapOriginal,obj.Smooth);
%             figure; tiledlayout('flow'); 
%             nexttile; imagesc(obj.MapOriginal);
%             nexttile; imagesc(obj.MapSmooth);
            obj.XEdges=Xedges;
            obj.ZEdges=Zedges;
            obj.Empties = [];
        end


        function obj = calculateAJ(obj,xedges,zedges)
            %adam version

            pd=obj.PositionData.data;
            sr = obj.PositionData.time.getSampleRate;
            x = pd.X; z = pd.Z;
            %Pos1=[pd.X pd.Z];
            x(belowspeedThresh) = 0;
            z(belowspeedThresh) = 0;
            binsize = obj.SpatialBinSizeCm;
            if nargin==1 % no xedges/zedges provided
                xbins = min(x):binsize:max(x);      % The x-coordinate bins for the grid
                zbins = min(z):binsize:max(z);      % The y-coordinate bins for the grid
            else
                xbins = xedges(1):binsize:xedges(2);
                zbins = zedges(1):binsize:zedges(2);
            end
            xbins = xbins(1:end-1) + 0.5*binsize;
            zbins = zbins(1:end-1) + 0.5*binsize;
            nxbins = length(xbins);
            nzbins = length(zbins);
            xbinInd = 1:nxbins;
            zbinInd = 1:nzbins;

            % Find the occupancy - how long does the animal spend in each location?
            occ = zeros(nxbins,nzbins);

            xb = interp1(xbins,xbinInd,x,'nearest');
            zb = interp1(zbins,zbinInd,z,'nearest');
            t = obj.PositionData.time.getTimePointsInSamples/sr; %in seconds
            xzt = [x z t'];
            speedthresh = 10; %cm per second. %%% YOU MIGHT WANT TO CHANGE THIS
            [~,belowspeedThresh] = external.velocityfilterUKRW_2D(xzt,speedthresh); %in cm/s
            
           
            for iX = 1:nxbins
                % Find all X positions in the bin
                xtmp = (xb == iX);
                for iZ = 1:nzbins
                    % Find all Y positions in the bin
                    ztmp = (zb == iZ);
                    occ(iX,iZ) = sr*sum(xtmp.*ztmp);
                end
            end
            empties = find(occ <= sr*10); %2 seconds minimum occupancy
      
            
            rr = 1;
            [xx,zz] = meshgrid(-2*rr:2*rr,-2*rr:2*rr);
            ff = exp(-((xx).^2 + (zz).^2)/4); %the filter
            ff = ff/(sum(ff(:)));
            % figure; surfl(xx,yy,zz)
            smoothedocc = filter2(ff,occ,'same');

           % occ(empties) = nan;
           % smoothedocc(empties) = nan;

            obj.MapOriginal=occ;
            obj.MapSmooth=smoothedocc;
            obj.XEdges=xbins;
            obj.ZEdges=zbins;
            obj.Empties=empties;
            
            %   figure;
            %     h = pcolor(smoothedocc');
            %     axis equal
            %     colorbar
            %     title('Time spent in each position')
            %     set(h,'linestyle','none')
            %     axis equal
  
        end

        function frm = plus(obj,spikeTimes)
            frm=neuro.placeField.FireRateMap(obj,spikeTimes);
        end
        function xl = getXLim(obj)
            xl=[min(obj.XEdges) max(obj.XEdges)];
        end
        function yl = getZLim(obj)
            yl=[min(obj.ZEdges) max(obj.ZEdges)];
        end
    end
end

