classdef MethodsFigure < FigurePipeline

    % a) example track
    % b) example traces
    % c) example histology electrode placement
    % d) raw example trace
    % e) th example trace
    % f) ripple example trace
    % g) ex raster
    % h) example place cells
    % i) violin plots for summary stats
    % j) more summary stats (m vs f?)
    % k) exp timeline - to be added in illustrator
    % Load behavioral data

    properties
        LFP
    end

    methods

        function obj = MethodsFigure()
            % to call, mf = MethodsFigure()

            obj.initializeProperties(); % also cd's to basepath

            %%%%%%%%%%%%%%%%% lfp %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            ctdh = neuro.basic.ChannelTimeDataHard(obj.basepath);
            lfpchannel = 30;
            obj.LFP = ctdh.getChannel(lfpchannel);
            % set zt time as noon
            ticd = obj.LFP.TimeIntervalCombined;
            ticd = ticd.setZeitgeberTime(hours(12)); % for noon start, dark to light transition
            obj.LFP.TimeIntervalCombined = ticd; %pulling the 1250Hz version!

        end

        function [] = exampleTrace(obj,timeframe,freqs)
            %time frame in seconds [start end]
            %freqs in hertz [low high]

            % manually selecting 60 minutes in to start (puts mid on track)
            % lfpStart = 60*60*1250;
            % lfpEnd = lfpStart + 10*1250; % 10 sec later
            sr = obj.TimeIntervalCombined.getSampleRate;
            lfpSeg = obj.LFP.Values(timeframe(1)*sr:timeframe(2)*sr);

            bandpassed = bandpass(lfpSeg,[5 12],sr);

            ax = gca;
            plot(bandpassed)
            %ax.XTick()
        end

        function ripples = getRipples(obj)
            %ripchan = bz_GetBestRippleChan_RW(cleanlfp); % most power in 140 to 180 Hz band. LFP has to come in as 1250 Hz
            ripCh = 126;
            noiseCh = 5;

            % Params (defaults listed below from buzcode tutorial
            % restrict periods chosen by plotting and selecting quiet period
            ripdur = [30 300];% was [30 200]  %100 was default, 400 KD's suggestion
            ripfreq = [140 240]; %ripfreq = [130 200]; %200 was default, 250 KD's suggestion

            premaze_ms = length(premazelfp)/1250*1000;
            maze_ms = length(mazelfp)/1250*1000;

            wpremazelfp = WhitenSignal(lfp.data(1:lfpstartIndex-1,ripCh));
            noiselfp = WhitenSignal(lfp.data(1:lfpstartIndex-1,noiseCh));
            timestamps = lfp.timestamps(1:lfpstartIndex-1);
            premazeRipples = bz_FindRipples(wpremazelfp,timestamps,'thresholds',[1 4],...
                'durations',ripdur,'show','on','saveMat',false,'noise',noiselfp,'EMGThresh',[]);
            clear wpremazelfp noiselfp tsindex timestamps
            disp('Finished finding ripples for premaze period')

            wmazelfp = WhitenSignal(lfp.data(lfpstartIndex:lfpendIndex,ripCh));
            noiselfp = WhitenSignal(lfp.data(lfpstartIndex:lfpendIndex,noiseCh));
            timestamps = lfp.timestamps(lfpstartIndex:lfpendIndex);
            mazeRipples = bz_FindRipples(wmazelfp,timestamps,'thresholds',[2 5],...
                'durations',ripdur,'show','on','saveMat',false,'noise',noiselfp,'EMGThresh',[], ...
                'restrict',[9000 max(timestamps)]); % identify wh
            clear wmazelfp noiselfp tsindex timestamps
            disp('Finished finding ripples for maze period')

            wpostmazelfp = WhitenSignal(lfp.data(lfpendIndex+1:end,ripCh));
            noiselfp = WhitenSignal(lfp.data(lfpendIndex+1:end,noiseCh));
            timestamps = lfp.timestamps(lfpendIndex+1:end);
            postmazeRipples = bz_FindRipples(wpostmazelfp,timestamps,'thresholds',[1 4],...
                'durations',ripdur,'show','on','saveMat',false,'noise',noiselfp,'EMGThresh',[], ...
                'restrict',[14000 18000]);
            clear wpostmazelfp noiselfp tsindex timestamps
            disp('Finished finding ripples for postmaze period')

            ripples = vertcat(...
                [premazeRipples.timestamps(:,1) premazeRipples.peaks(:,1) premazeRipples.timestamps(:,2)], ...
                [mazeRipples.timestamps(:,1) mazeRipples.peaks(:,1) mazeRipples.timestamps(:,2)], ...
                [postmazeRipples.timestamps(:,1) postmazeRipples.peaks(:,1) postmazeRipples.timestamps(:,2)]);

            SaveRippleEvents_RW(fullfile(basepath,[basename '.rip.evt']),ripples,ripCh) %ripchan is 127

            ripples = ripples*1000; % to convert into ms

        end

        function [] = testAlignment(obj)

            % test alignment of ZT timestamps across spikes, position

            % pull data for the specified time window (see figureScript.m)
            timeWindow = [obj.positionData.time.getStartTimeAbs + minutes(45) obj.positionData.time.getStartTimeAbs + minutes(48)];
            mftw = MethodsFigureTimeWindow(obj,timeWindow); %otherwise MethodsFigure gets overwritten!

            ztPosition = mftw.positionData.time.getTimePointsZT;
            [powerTheta,tfmTheta] = mftw.getPower([5 12]);
            [powerSWR,tfmSWR] = mftw.getPower([150 250]);
            s = mftw.positionData.getSpeed; %
            svals = s.Values;
            sbad = find(svals>200);
            svals(sbad) = nan;

            fullWindow(1) = obj.LFP.TimeIntervalCombined.getStartTimeAbs;
            ti4 = obj.LFP.TimeIntervalCombined.timeIntervalList.get(4);
            fullWindow(2) = ti4.getEndTimeAbs;
            obj.getStates(fullWindow);
            % obj.getStates(timeWindow);

            figure
            hold all
            grid_height = 6; grid_width = 1;
            h1 = tiledlayout(grid_height,grid_width);
            position = 1; h = 1; w = 1; ax1 = nexttile(position,[h,w]); %position
            position = 2; h = 2; w = 1; ax2 = nexttile(position,[h,w]); %raster
            position = 4; h = 1; w = 1; ax3 = nexttile(position,[h,w]); %
            position = 5; h = 1; w = 1; ax4 = nexttile(position,[h,w]); % acg whole session
            position = 6; h = 1; w = 1; ax5 = nexttile(position,[h,w]);
            % position = 6; h = 1; w = 1; ax6 = nexttile(position,[h,w]);
            %position = 7; h = 1; w = 1; ax7 = nexttile(position,[h,w]);

            axes(ax1); % X Position
            hold on
            xlabel('Time')
            ylabel('Position')
            plot(ztPosition,mftw.positionData.data.X,'r')
            plot(ztPosition,mftw.positionData.data.Z,'b')
            % plot(ztPosition,mftw.positionData.data.Y,'g')
            %add a legend here too
            hold off

            axes(ax2); % spike raster (just for task epoch)
            hold on
            xlabel('Time')
            ylabel('Cell')
            mftw.spikeRaster(timeWindow);
            % ylim([ztPosition(1) ztPosition(end)])
            hold off

            axes(ax3); % theta power
            hold on
            xlabel('Time')
            ylabel('Frequency')
            tfmTheta.plot('ZT');
            % ylim([timeWindow(1) timeWindow(2)]);

            axes(ax4); %SWR power
            hold on
            xlabel('Time')
            ylabel('Frequency')
            tfmSWR.plot('ZT')

            axes(ax4)
            hold
            xlabel('Time')
            ylabel('uV')
            plot(swrband.TimeIntervalCombined.getTimePointsZT,swrband.Values)

            axes(ax5);
            hold on
            xlabel('Time')
            ylabel('Speed (cm/s)')
            plot(ztPosition,svals)
            ylim([min(svals) max(svals)])

            linkaxes([ax1 ax2 ax3 ax4 ax5],'x')
            %  ax1.Visible="off";
            % ax2.Visible="off";
            %ax3.Visible="off";

            h1.TileSpacing='tight';
        end

        function eegStates = getStates(obj,timeWindow)

            %%%%%%%%%%%%%Check EEG states with sliding window %%%%%%%%%%%%%%%%%%%%%%%%%
            % Put "100" as a window size into WhitenSignal, didn't change much
            % figure out what the measurement is for the window

            %%%%%%%%%%%%%Check EEG states for entire behavioral period %%%%%%%%%%%%%%%%
            lfp = obj.LFP.getTimeWindow(timeWindow);

            eegStates = neuro.power.CheckEegStates(obj.basepath,obj.basename,'lfp',lfp.Values','redo_flag',true);
            %select ch 1 because I'm just sending in 1 channel anyway

            %%%%%%%%%%%%%% Try SleepScoreMaster for theta detection %%%%%%%%%%%%%%%%%%

            %badCh = session.channelTags.Bad.channels;
            %thetaCh = session.channelTags.ThetaChan.channels; % check this

            % SleepState = SleepScoreMaster(basepath,'ignoretime',removetimestamps.sampformat,'rejectchannels',badCh,'ThetaChannels',thetaCh);
        end

        function [pow,tfm] = getThetaPower(obj,timeWindow)

            % ctdh=ses.getDataLFP.getChannelTimeDataHard;
            % thetaLFP=ctdh.getChannel( ...
            %     ses.getDataLFP.getStateDetectionData.SleepScoreLFP.THchanID);
            thetaLFP = obj.LFP;
            thetaLFPt=thetaLFP.getTimeWindow(timeWindow);%(ses.getBlock('TRACK'));

            %             ticd = thetaLFPt.TimeIntervalCombined;
            %             ticd = ticd.setZeitgeberTime(hours(12)); % for noon start, dark to light transition
            %             thetaLFPt.TimeIntervalCombined = ticd; %pulling the 1250Hz version!

            tfm=thetaLFPt.getTimeFrequencyMap(...
                neuro.timefrequency.TimeFrequencyWavelet(5:.5:10));
            % ticd=time.TimeIntervalCombined(fullfile(theFile.folder, theFile.name)); % from the time folder
            pow=tfm.getPower(5:.5:10);
        end

        function [pow,tfm] = getPower(obj,bandpass)
            %timeWindow is easiest in absolute time, I think you can also
            %put it in as samples or zt but I was running into issues.
            %bandpass for the frequency band [lo hi] you would like.

            lfp = obj.LFP.getBandpassFiltered(bandpass);

            % set zt time as noon
            %             ticd = lfp.TimeIntervalCombined;
            %             ticd = ticd.setZeitgeberTime(hours(12)); % for noon start, dark to light transition
            %             lfp.TimeIntervalCombined = ticd; %pulling the 1250Hz version!

            % get time frequency map
            tfm=lfp.getTimeFrequencyMap(...
                neuro.timefrequency.TimeFrequencyWavelet(bandpass(1):.5:bandpass(2)));
            % ticd=time.TimeIntervalCombined(fullfile(theFile.folder, theFile.name)); % from the time folder
            pow=tfm.getPower(bandpass(1):.5:bandpass(2));

        end

        function [pow,tfm] = getPowerTimeWindow(obj,timeWindow,bandpass)
            %timeWindow is easiest in absolute time, I think you can also
            %put it in as samples or zt but I was running into issues.
            %bandpass for the frequency band [lo hi] you would like.
            lfp = obj.LFP.getTimeWindow(timeWindow);
            lfp = lfp.getBandpassFiltered(bandpass);

            %set zt time
            ticd = lfp.TimeIntervalCombined;
            disp('input zt for getPowerTimeWindow')
            prompt = {'Zeitgeber Time:'};
            dlgtitle = 'title';
            dims = [1 10];
            definput = {'12:00'};
            zt = duration(inputdlg(prompt,dlgtitle,dims,definput),'InputFormat','hh:mm');
            ticd = ticd.setZeitgeberTime(zt); % for noon start, dark to light transition
            lfp.TimeIntervalCombined = ticd; %pulling the 1250Hz version!

            % get time frequency map
            tfm=lfp.getTimeFrequencyMap(...
                neuro.timefrequency.TimeFrequencyWavelet(bandpass(1):.5:bandpass(2)));
            % ticd=time.TimeIntervalCombined(fullfile(theFile.folder, theFile.name)); % from the time folder
            pow=tfm.getPower(bandpass(1):.5:bandpass(2));
        end

        function [] = plotSpeedandRipples



        end

        function [] = spikeRaster(obj,timeWindow)

            %%first argument is what column to sort into different figures, second arg
            %is the filepath to which to save figures (kills matlab with too many
            %units!)

            sf = neuro.spike.SpikeFactory.instance(); %Spike Factory (Kaya code)
            sa = sf.getPhyOutputFolder(obj.phypath); %Spike Array (Kaya code)

            if exist("timeWindow","var") % pull out timeWindow, otherwise takes the whole thing
                sa = sa.getTimeWindow(timeWindow);
            end

            sorted = sa.sort('fr');

            sav = neuro.spike.SpikeArrayVisualize(sorted); %class with visualization functions
            % is currently creating a figure of its own
            sav.plotRasterZT('class',[1]); %('sh','/home/wahlberg/merged_20230614crs.GUI/Figures/')
        end

        function plotPlaceCellsKD(obj)

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

            %obj.positionData = obj.positionData.normalizePositionData;
            %%%%%%%%% filter by... PART 1 %%%%%%%%%%%

            %%% Clockwise + Counterclockwise indices
            %[clockwiseIdx,counterIdx] = getDirectionsIdx(obj.Events.Turns);
            %positionTimeZT = positionTimeZT(clockwiseIdx);

            %%% First half/ second half 
            half = 'Both'; %'First' 'Second'
            %[positionTimeZT,normX,normZ] = filterbyhalf(positionTimeZT,normX,normZ,half);

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

            %sat = neuro.spike.SpikeArrayTrack(sa,obj.positionData); %SpikeArrayTrack (Kaya)
            %sut = obj.getPlaceFields(sat); %get just the spike info for the track period

            %%%%%% for each unit, plot place field %%%%%%%
            figure
            grid_height = 2; grid_width = 3;
            for c = 1:nclassTypes
                class = classTypes(c);
                for unit = 1:nUnits

                    %%%%%% get spikeCount %%%%%%%%
                    % spsBinned is binned with same bins as binnedPos above

                    spksZT = seconds(su(unit).getTimesZT)';  %spikes in ZT time.

                    spksBinned = histcounts(spksZT,histedges)';
                    spksBinned(end)=0; spksBinned(1)=0; % find number of spikes per unit time.

                    spkslogical = logical(spksBinned);
                    xprevelocityfilter = positionBinned(spkslogical,1)*100;
                    zprevelocityfilter = positionBinned(spkslogical,2)*100;

                    %%%%%%%%% % velocityfilter = external.velocityfilter(binnedPos)
                    velocityThresh = 10; %in 10 cm/hour
                    [~,~,spksBinned] = external.velocityfilterUKRW_2D(positionBinned,velocityThresh,spksBinned); %in cm/s
                    
                    % spksBinned(belowthreshTheta) = 0;
                    % Set limits defined by the shape of the track
                    % limits.x = [0.3 0.7];
                    % limits.z = [0.3 0.7];

                    %  [xztRange, nspkRange] = position.SetCircularRange(binnedPos,spkCount2,limits);

                    spkslogical = logical(spksBinned);
                    xpostvelocityfilter = positionBinned(spkslogical,1)*100;
                    zpostvelocityfilter = positionBinned(spkslogical,2)*100;

                    h1 = tiledlayout(grid_height,grid_width);
                    title(h1,['Unit ' num2str(unit)]);
                    position = 1; h = 1; w = 1; ax1 = nexttile(position,[h,w]);
                    position = 2; h = 1; w = 1; ax2 = nexttile(position,[h,w]);
                    position = 3; h = 1; w = 1; ax3 = nexttile(position,[h,w]); 
                    position = 4; h = 1; w = 1; ax4 = nexttile(position,[h,w]);
                    position = 5; h = 1; w = 1; ax5 = nexttile(position,[h,w]);
                    position = 6; h = 1; w = 1; ax6 = nexttile(position,[h,w]); 

                    % AX1: plots position data vertically over time, and then the locations of unit firing
                    axes(ax1)
                    plottitle = 'Spk locs. (all velocities)';
                    plot.plotPlaceFields('spkLocations',plottitle,x,z,xprevelocityfilter,zprevelocityfilter)

                    % AX2: plots position data vertically over time, and then the locations of unit firing
                    axes(ax2)
                    plottitle = 'Spk locs. (> 10 cm velocities)';
                    plot.plotPlaceFields('spkLocations',plottitle,x,z,xpostvelocityfilter,zpostvelocityfilter)

                    % AX3: Plots occupancy and firing rate map
                    axes(ax6)
                    hold on
                    title(['Place map, nspikes = ' num2str(sum(sum(nSpikes)))])
                    %Get and plot place field %%%%%%%
                    [~,~,nSpikes,sTimeSpent,snSpikes] = ...
                        external.Placefields.PFClassic(positionBinned(:,1:2),spksBinned,smoothVal,nGrid,timebinsize);

                    % AX4:
                    axes(ax4)
                    plottitle = 'Smoothed Occ. map';
                    plot.plotPlaceFields('imagescMap',plottitle,sTimeSpent)

                    % AX5:
                    axes(ax5)
                    plottitle = 'Smoothed Spikemap';
                    plot.plotPlaceFields('imagescMap',plottitle,snSpikes)

                    pause
                    clf
                end
            end
        end

        function [] = checkOccupancy(obj)
            % refer to pc_clockcounter_linear or similar
            import position.* %otherwise when I call position.getDirectionsIdx it's calling a different position func?

            %%%%%% default vals %%%%%%%%
            %timeWindow = [mf.positionData.time.getStartTimeAbs mf.positionData.time.getEndTimeAbs];
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

            %Get binned position %%%%%%%%%%%
            tbegin = positionTimeZT(1); tend = positionTimeZT(end);
            nBins = round(seconds(tend-tbegin)/timebinsize);
            xzt(:,1) = normX;
            xzt(:,2) = normZ;
            xzt(:,3) = seconds(positionTimeZT);

            positionBinned = external.Placefields.binpos(xzt,nBins);
            histedges = positionBinned(:,3)-timebinsize;
            histedges(end+1) = histedges(end)+timebinsize;

            %%%%%%%%% % velocityfilter = external.velocityfilter(binnedPos)
            velocityThresh = 10; %in 10 cm/hour
            goodindex = velocityfilter_2D(positionBinned,velocityThresh); %in cm/s (in +position folder)

            %%%% get angle
            [Turns,angle] = directionsCircularTrackRadians(obj.positionData);
            at = [angle xzt(:,3)];
            angleBinned=external.Placefields.binAngle(at,nBins); %binned at same rate as binnedPosition

             allBins = ones(length(positionBinned),1);
            allAngles = getSpksAngle(allBins,positionBinned);
  
            %%%%% check occupancy %%%%
            abovethreshAll = PlaceMaps.peakAll.FiringRate > frThresh;
            abovethreshPlaceMaps = PlaceMaps.placeMapAll(abovethreshAll,:);
            abovethreshTimeSpent = PlaceMaps.timeSpentAll(abovethreshAll,:);
            abovethreshPeaks = PlaceMaps.peakAll(abovethreshAll,:);

            [~,positionOrder] = sort(abovethreshPeaks.Position,'ascend');
            abovethreshPlaceMaps = abovethreshPlaceMaps(positionOrder,:);
            abovethreshTimeSpent = abovethreshTimeSpent(positionOrder,:);

            neuro.placeField.PFPlot_RW(abovethreshPlaceMaps,abovethreshTimeSpent, TopRate, TimeThresh);
       
            
            allAngles(allAngles<0) = allAngles(allAngles<0)+(2*pi);%to convert from -pi to pi range to 0 to 2pi 
            %allAngles = allAngles/max(allAngles); %not subtracting min first on purpose, because I want to keep true zero

            edges = linspace(0,2*pi,101);
            [N,edges] = histcounts(allAngles,edges);
            figure
            histogram('BinEdges',edges,'BinCounts',N)

            aboveVelocityThresh = allAngles;
            aboveVelocityThresh(~goodindex) = 10; %out of range of 0 to pi
            [N,edges] = histcounts(aboveVelocityThresh,edges);
            figure
            histogram('BinEdges',edges,'BinCounts',N)
            xlabel('Radians')
            ylabel('Bin Count')
            title('Linearized Occupancy - above 10 cm/s')



        end

    end

    %%%%%  OTHER PLOT COMBINATIONS

    methods static %these are the SAME CODE, just filtered in different ways
        function [] = plotPlaceCellsKD_velocitycomparison(obj) %all velocities vs >10cm only
            neuro.placeField.pc_velocityfiltercomparison;
        end

        function [] = plotPlaceCellsKD_clockwiseCounter(obj) %compare counter vs clockwise plots
            neuro.placeField.pc_clockcountercomparison;
        end

        function [] = plotPlaceCellsKD_clockwiseCounterLinear(obj) %compare counter vs clockwise plots
            plotfigures = 1;
            PlaceMaps = neuro.placeField.pc_clockcounter_linear(obj,plotfigures);
          
            frThresh = 1;
            TopRate = [];
            TimeThresh = 0.2;

            abovethreshAll = PlaceMaps.peakAll.FiringRate > frThresh;
            abovethreshClustersAll = obj.ClusterInfo(abovethreshAll,:);
            belowthreshClustersAll = obj.ClusterInfo(~abovethreshAll,:);
            abovethreshPlaceMaps = PlaceMaps.placeMapAll(abovethreshAll,:);
            abovethreshTimeSpent = PlaceMaps.timeSpentAll(abovethreshAll,:);
            abovethreshPeaks = PlaceMaps.peakAll(abovethreshAll,:);

            [~,positionOrder] = sort(abovethreshPeaks.Position,'ascend');
            abovethreshPlaceMaps = abovethreshPlaceMaps(positionOrder,:);
            abovethreshTimeSpent = abovethreshTimeSpent(positionOrder,:);
            
            neuro.placeField.PFPlot_RW(abovethreshPlaceMaps,abovethreshTimeSpent, TopRate, TimeThresh);
            
            %%%%% clockwise
            abovethreshClock = PlaceMaps.peakClock.FiringRate > frThresh;
            abovethreshClustersClock = obj.ClusterInfo(abovethreshClock,:);
            belowthreshClustersClock = obj.ClusterInfo(~abovethreshClock,:);

            abovethreshPlaceMaps = PlaceMaps.placeMapClock(abovethreshClock,:);
            abovethreshTimeSpent = PlaceMaps.timeSpentClock(abovethreshClock,:);
            abovethreshPeaks = PlaceMaps.peakClock(abovethreshClock,:);

            [~,positionOrder] = sort(abovethreshPeaks.Position,'ascend');

            abovethreshPlaceMaps = abovethreshPlaceMaps(positionOrder,:);
            abovethreshTimeSpent = abovethreshTimeSpent(positionOrder,:);
            
            neuro.placeField.PFPlot_RW(abovethreshPlaceMaps,abovethreshTimeSpent, TopRate, TimeThresh);

            %%%% counter

            abovethreshCounter = PlaceMaps.peakCounter.FiringRate > frThresh;
            abovethreshClustersCounter = obj.ClusterInfo(abovethreshCounter,:);
            belowthreshClustersCounter = obj.ClusterInfo(~abovethreshCounter,:);

            abovethreshPlaceMaps = PlaceMaps.placeMapCounter(abovethreshCounter,:);
            abovethreshTimeSpent = PlaceMaps.timeSpentCounter(abovethreshCounter,:);
            abovethreshPeaks = PlaceMaps.peakCounter(abovethreshCounter,:);

            [~,positionOrder] = sort(abovethreshPeaks.Position,'ascend');

            abovethreshPlaceMaps = abovethreshPlaceMaps(positionOrder,:);
            abovethreshTimeSpent = abovethreshTimeSpent(positionOrder,:);
            
            neuro.placeField.PFPlot_RW(abovethreshPlaceMaps,abovethreshTimeSpent, TopRate, TimeThresh);

            figure
            plot(obj.positionData.data.X,obj.positionData.data.Z,'.')
            hold on
            yline(0.5)
            xline(0.5)
            line([0,1],[0,1])
            line([0,1],[1,0])

            % 3/11  @ 0 - 0
            % 2/10  @ pi/4 - 0.785
            % 1/9  @ pi/2 -  1.571
            %  0/8  @ 3pi/4 - 2.36
            %  7/15 @ pi - 3.14
            %  6/14@ 5pi/4 - 3.93
            %  5/13 @ 3pi/2 - 4.71
            % 4/12 @ 7pi/4 - 5.45
        end

    end

    methods

        function [sut,acgTask,frm,pfm] = getPlaceFields(obj,sa)
            spikeUnits = sa.getSpikeUnits();

            try
                classTypes = unique(sa.ClusterInfo.class);
            catch
                classTypes = 1; %just for UK data
            end

            for c = 1:length(classTypes)
                class = classTypes(c);
                unitcount = 1;
                for unit = 1:length(spikeUnits)
                    if spikeUnits(unit).Info.class == class
                        try
                            sut{c}{unitcount} = neuro.spike.SpikeUnitTracked(spikeUnits(unit),obj.positionData);
                            acgTask{c}{unitcount} = sut{c}{unitcount}.getACC;
                            frm{c}{unitcount} = sut{class}{unitcount}.getFireRateMap; %underlay is occupancy map
                            pfm{c}{unitcount} = frm{c}{unitcount}.getPlaceFieldMap;
                        catch
                            sut{c}{unitcount} = [];
                            acgTask{c}{unitcount} = [];
                            frm{c}{unitcount} = [];
                            pfm{c}{unitcount} = [];
                        end
                        unitcount = unitcount + 1;
                    end % if spikeUnits
                    unit
                end % for unit =
            end %for c =
        end %function

        function [sut,acg] = getEntireSession(obj,sa)

            % get ACGs and waveforms BEFORE grabbing out just the task ep
            spikeUnits = sa.getSpikeUnits();
            classTypes = unique(sa.ClusterInfo.class);
            for c = 1:3%length(classTypes)
                class = classTypes(c);
                unitcount = 1;
                for unit = 1:length(spikeUnits)
                    if spikeUnits(unit).Info.class == class
                        try
                            sut{c}{unit} = neuro.spike.SpikeUnitTracked(spikeUnits(unit),obj.positionData);
                            acg{c}{unitcount} = sut{c}{unit}.getACC;
                        catch
                            sut{c}{unit} = [];
                            acg{c}{unitcount} = [];
                        end
                        unitcount = unitcount + 1;
                    end
                    unit
                end
            end
        end

        function [] = testSubplots(obj,toplot) % so you can test ind. figures,
            %or plot into a full figure with plotfullfig

            figure

            if strcmp(toplot,'exampleTraces') == 1
                exampleTraces(obj)
            end

            if strcmp(toplot,'spikeRaster') == 1
                spikeRaster(obj)
            end

        end

        function [] = summaryStats(obj)

        end


        function [] = plotFigure(obj)

            positionData = obj.loadPositionData;

            h1 = figure;
            grid_height = 6; grid_width = 6;
            tiledlayout(grid_height,grid_width);

            %%a - example track: add in illustrator
            % position = 1; height = 2; width = 2; nexttile(position,[height,width]);

            %%b example behavior
            position = 13; height = 1; width = 1; ax1 = nexttile(position,[height,width]);


            plot(positionData.data.X,positionData.data.Z,'LineWidth',2)

            %%c - histology placement: add in illustrator
            % position = 14; height = 1; width = 1; nexttile(position,[height,width]);

            %%d example traces
            position = 3; height = 3; width = 4; ax2 = nexttile(position,[height,width]);
            tiledlayout(height,width);
            exampleTraces(obj)

            %%e theta example trace
            %position = 9; height = 1; width = 4; nexttile(position,[height,width]);

            %%f ripple example trace
            %position = 15; height = 1; width = 4; nexttile(position,[height,width]);

            %%g example raster
            position = 21; height = 1; width = 4; ax3 = nexttile(position,[height,width]);
            % spikeRaster(obj)

            %%h example placecells (includes nested tiledlayout)
            position = 19; height = 2; width = 2; ax4 = nexttile(position,[height,width]);

            % neuro.spike.spikeUnit

            %             spikeIDs = sa.getspikeIDs;
            %             st = sa.getSpikeTimes;
            %             ticd = neuro.time.TimeIntervalCombined(phyfolder);
            %
            %             su = neuro.spike.SpikeUnit(spikeIDs,st,ticd); %


            %%i violin plots for summary stats
            position = 27; height = 1; width = 2; ax5 = nexttile(position,[height,width]);

            %%j more summary stats?
            position = 29; height = 1; width = 2; ax6 = nexttile(position,[height,width]);

        end

        function [] = saveFigure(filename)
            plotFigure()
            saveas(gcf,filename)
        end
    end

end






