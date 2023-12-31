classdef PlaceCells < FigurePipeline
    % for place cell analyses

    properties
        fBase = 1.0;
        minT_occ = 1.5;
        tau = 0.25;
        binsize = 2; % cm
        %%%% inherited properties %%%%
%         basepath = '/data/ExperimentsRSW/CircularMaze/Harry/CA1/20230614/merged_20230614';
%         basename = 'merged_20230614';
%         behaviorpath = '/data/ExperimentsRSW/CircularMaze/Harry/CA1/20230614/merged_20230614/behaviorfiles/';
%         phypath = '/home/wahlberg/merged_20230614crs.GUI/';
%         optipath = '/data/ExperimentsRSW/CircularMaze/Harry/CA1/20230614/Optitrack';
%         SpikeTable
%         ClusterInfo
%         Info
%         positionData
%         TimeIntervalCombined
    end

    methods (Static)

        function obj = PlaceCells(obj)
            obj.initializeProperties();
        end

        function [] = plot2Dfield_oneunit(unitID,timewindow)
            %time window - can select small windows to observe time change or entire epoch, for ex

            
            lfpSR = obj.TimeInterValCombined.getSampleRate;
            behSR = obj.BehaviorTICD.getSampleRate;

            interpolatedData = obj.interpPositionData;
            x = interpolatedData.X(timewindow(1)*lfpSR:timewindow(2)*lfpSR);
            y = interpolatedData.Y(timewindow(1)*lfpSR:timewindow(2)*lfpSR);

            xbins = min(x):binsize:max(x);      % The x-coordinate bins for the grid
            ybins = min(y):binsize:max(y);      % The y-coordinate bins for the grid
            xbins = xbins(1:end-1) + 0.5*binsize;
            ybins = ybins(1:end-1) + 0.5*binsize;
            nxbins = length(xbins);
            nybins = length(ybins);
            xbinInd = 1:nxbins;
            ybinInd = 1:nybins;

            % Find the occupancy - how long does the animal spend in each location?
            occ = zeros(nxbins,nybins);
            xb = interp1(xbins,xbinInd,x,'nearest');
            yb = interp1(ybins,ybinInd,y,'nearest');
            for iX = 1:nxbins
                % Find all X positions in the bin
                xtmp = (xb == iX);
                for iY = 1:nybins
                    % Find all Y positions in the bin
                    ytmp = (yb == iY);
                    occ(iX,iY) = behsamprate*sum(xtmp.*ytmp);
                end
            end
            empties = find(occ < 1);

            %%%% Smoothing kernel
            %   We blur the position of the animal because we don't exactly where it
            %   is.  The smoothing kernal basically spreads the position of the animal
            %   at each video tracker time as a 2D Gaussian (normal) distribution with
            %   a center where we believe the animal's head was according to the video.
            %   4.2 pixel/cm
            rr = 1;
            [xx,yy] = meshgrid(-2*rr:2*rr,-2*rr:2*rr);
            zz = exp(-((xx).^2 + (yy).^2)/4); %the filter
            zz = zz/(sum(zz(:)));
            % figure; surfl(xx,yy,zz)
            smoothedocc = filter2(zz,occ,'same');
            smoothedocc(empties) = nan;

            %             % prior, probability pre-data of being in the group
            % posterior, probability of being in the group given the data
            %
            alpha0 = fBase*minT_occ;    % alpha hyperparameter for gamma prior
            beta0 = minT_occ/tau;       % beta hyperparameter for gamma prior
            alpha = alpha0 + smoothedspike;   % alpha hyperparameter for gamma posterior
            beta = beta0 + smoothedocc/tau;   % beta hyperparameter for gamma posterior
            likelihood = alpha; % probability of spikes
            marg = beta./(beta + 1); % marginilization % prob of occupancy
            posterior = likelihood.*(1-marg)./marg; % p(A|B) = p(B|A).*(P(A)./p(B)    % called meanN
            %prob of alpha given beta, or prob of spikes given occupancy
            varN = likelihood.*(1-marg)./marg.^2;

            FRmap = smoothedspike./smoothedocc;     % This is your standard spikes per unit time based map.  It's not statistically correct.  But it's close and easy.
            FRmap = FRmap';

            %%% plot
            gcf % have a figure already in existence
            grid_height = 2; grid_width = 2;
            tiledlayout(grid_height,grid_width);

            nexttile(1,[1 1]);
            h = pcolor(FRmap);
            set(h,'linestyle','none')
            axis equal
            %         caxis([0 40]);
            colorbar
            axis off
            m = find(FRmap == max(FRmap(:)));
            [xtmp ytmp] = ind2sub(size(FRmap),m);
            hold on;
            plot(ytmp,xtmp,'wo','markersize',15);
            hold off;
            title('Raw mean FR (Hz)','fontsize',15);

            % Make spiking plot
            nexttile(2,[1 1]);
            plot(x,y,'color',[0.5 0.5 0.5]);
            hold on;
            plot(x,y,'r.');
            title(['Unit: ' num2str(iU)],'fontsize',15);
            colorbar
            axis equal
            axis off

            nexttile(3,[1 1]);
            hold on
            plot(t,xb);
            if ~isempty(x); scatter(1:length(x),x,'r','LineWidth',2); end;
            ylabel('X pos')
            xlabel('time')
            title('X versus time')

            nexttile(4,[1 1]);
            hold on
            plot(t,yb);
            if ~isempty(y); scatter(1:length(y),y,'r','LineWidth',2); end;
            ylabel('Y pos')
            xlabel('time')
            title('Y versus time')
        end

        function [] = plotall2Dfields(group,id)
            %allows you to plot all pyramidals, for instance; pauses in between

            % sort by group,labels within the group


            % plot all
            figure
            for iU = 1:nUnits
                plot2Dfield_oneunit(iU)
                title(['Unit ' num2str(iU) ', ' type])
                pause
            end
        end






    end %methods end
end %classdef end


