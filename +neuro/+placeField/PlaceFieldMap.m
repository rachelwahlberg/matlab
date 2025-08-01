classdef PlaceFieldMap<neuro.placeField.FireRateMap
    %PLACEFIELDMAP Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        TimeThreshold
        FireRateMap
        Parent
    end

    methods
        function obj = PlaceFieldMap(fireRateMap)
            %takes firing rate map (spike occupancy/position occupancy)
            if nargin>0
                fnames=fieldnames(fireRateMap);
                for ifn=1:numel(fnames)
                    obj.(fnames{ifn})=fireRateMap.(fnames{ifn});
                end
                obj.FireRateMap=fireRateMap.MapSmooth;
              %  obj.FireRateMap(obj.OccupancyMap<seconds(obj.TimeThreshold))=nan;

                obj.TimeThreshold=1*fireRateMap.TimeBin;

                fBase = 1.0; % The assumed basal firing rate for each unit
                minT_occ = 1; % Minimum sampling time (in seconds) at each location
                tau = 0.25; % The time-bin size (in seconds)

                % prior, probability pre-data of being in the group
                % posterior, probability of being in the group given the data
                %
                alpha0 = fBase*minT_occ;    % alpha hyperparameter for gamma prior
                beta0 = minT_occ/tau;       % beta hyperparameter for gamma prior
                alpha = alpha0 + fireRateMap.SmoothedSpike;   % alpha hyperparameter for gamma posterior
                beta = beta0 + fireRateMap.OccupancyMapSmooth/tau;   % beta hyperparameter for gamma posterior
                likelihood = alpha; % probability of spikes
                marg = beta./(beta + 1); % marginilization % prob of occupancy
                posterior = likelihood.*(1-marg)./marg; % p(A|B) = p(B|A).*(P(A)./p(B)    % called meanN
                %prob of alpha given beta, or prob of spikes given occupancy
                varN = likelihood.*(1-marg)./marg.^2;
        
                obj.MapSmooth = posterior'/tau;
%                 subplot(2,2,4)
%                 h = pcolor(fireRateMap.MapSmooth' - (posterior'/tau));
%                        %  h = pcolor(sqrt(varN)'/tau);
%                 set(h,'linestyle','none')
%                 axis equal
%                 colorbar
%                 axis off
%                 title('Raw Mean FR - NB-based mean FR (Hz)','fontsize',15);
%                 colormap('jet')
% 
%                 obj.MapSmooth=obj.FireRateMap./(obj.OccupancyMap+eps);
%                 obj.MapSmooth=imgaussfilt(obj.MapSmooth,obj.Smooth);
             end
        end

        function str=toString(obj)
            str=obj.SpikeUnitTracked.tostring;
        end

            function [] = plot(obj)
                h = pcolor(obj.MapSmooth);
                set(h,'linestyle','none')
        %        axis equal
         %       colorbar
          %      axis off
                title('NB-based mean FR (Hz)');%,'fontsize',15

            end

        function [] = plotUK(obj)
            ms=obj.OccupancyMap;
            ms(ms<eps)=0;
            alpha1=log(ms);
            x=[min(obj.PositionData.data.X) max(obj.PositionData.data.X)];
            y=[min(obj.PositionData.data.Z) max(obj.PositionData.data.Z)];
            imagesc(x,y,obj.MapOriginal,AlphaDataMapping="scaled",AlphaData=alpha1);
            xlabel(['X ' obj.Units])
            ylabel(['Z ' obj.Units])
        end
        function [] = plotSmooth(obj)
            ms=obj.OccupancyMap;
            if min(size(ms))>1
                is2d=1;
            else
                is2d=0;
            end
            ms(ms<eps)=0;
            alpha1=log(ms);
            x=[min(obj.PositionData.data.X) max(obj.PositionData.data.X)];
            y=[min(obj.PositionData.data.Z) max(obj.PositionData.data.Z)];
            if is2d
                imagesc(x,y,obj.MapSmooth,AlphaDataMapping="scaled",AlphaData=alpha1)';
            else
                imagesc(x,y,obj.MapSmooth)';
            end
            
            %             ax=gca;ax.CLim=[.05 .2];
            xlabel(['X ' obj.Units])
            ylabel(['Z ' obj.Units])
            str=sprintf('Info: %.3f bits\n',obj.Information);
            if ~isempty(obj.Stability)
                str=sprintf('%sStability: %.3f\n',str,obj.Stability.gini);
            end
            text(.5,.1,str, ...
                Units="normalized", ...
                VerticalAlignment="middle",HorizontalAlignment="center");
        end
        function [peak] = getPeakLocalMaxima(obj)
            [pks,locs1,w,p]=findpeaks(obj.MapSmooth);
            X=linspace(min(obj.PositionData.data.X), ...
                max(obj.PositionData.data.X), ...
                numel(obj.MapSmooth))';
            locs=X(locs1);
            peak=table(pks',locs,w',p', ...
                VariableNames={'FiringRate','Position','Width','Prominence'});
            peak=sortrows(peak,{'FiringRate', 'Prominence'},"descend");
        end
        function [peak] = getPeak(obj)
            [pks,locs1]=max(obj.MapSmooth);
             Z=linspace(min(obj.PositionData.data.Z), ...
                max(obj.PositionData.data.Z), ...
                size(obj.MapSmooth,2));
            locs=Z(locs1);
%             X=linspace(min(obj.PositionData.data.X), ...
%                 max(obj.PositionData.data.X), ...
%                 size(obj.MapSmooth,2));
%             locs=X(locs1);
            peak=table(pks',locs', ...
                VariableNames={'FiringRate','Position'});
        end
        function [ret] = getPlaceFieldMapMeasures(obj)
            ret=neuro.placeField.PlaceFieldMapMeasures(obj);
        end
        function [cache,obj] = getPositionHeld(obj,cache)
            [cache,key]=cache.hold( ...
                obj.PositionData);
            obj.PositionData=key;
            try
                [cache, obj.Parent]=obj.Parent.getPositionHeld(cache);
            catch ME
                
            end
            try
                [cache, obj.SpikeUnitTracked]=...
                    obj.SpikeUnitTracked.getPositionHeld(cache);
            catch ME
                
            end
        end
        function [obj] = getPositionReload(obj,cache)
            try
                obj1=cache.get(obj.PositionData);
                obj.PositionData=obj1;
            catch ME
                
            end
            try
                obj.Parent=obj.Parent.getPositionReload(cache);
            catch ME
                
            end
            try
                obj.SpikeUnitTracked=obj.SpikeUnitTracked.getPositionReload(cache);
            catch ME
                
            end
        end
    end
end

