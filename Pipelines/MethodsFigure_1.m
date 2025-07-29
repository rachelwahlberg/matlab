classdef MethodsFigure_1 < FigurePipeline

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

        function obj = MethodsFigure_1()
            % to call, mf = MethodsFigure()
            
            obj.initializeProperties(); % also cd's to basepath

            %%%%%%%%%%%%%%%%% lfp %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            ctdh = neuro.basic.ChannelTimeDataHard(obj.basepath);
            lfpchannel = 30;
            obj.LFP = ctdh.getChannel(lfpchannel);
              % set zt time as noon
            ticd = obj.LFP.TimeIntervalCombined;
            ticd = ticd.setZeitgeberTime(hours(8)); % for noon start, dark to light transition
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

        function [] = plotPlaceMaps1(obj,idx) %Occupancy currently looks cut off!
            %not fully working Dec 2023 RW

            % each place cell in a cell
            %idx = 'class';
            import neuro.spike.*
            import neuro.placeField.*

            % pd_ontrack = obj.positionData.getTimeWindow(timewindow)
            sf = neuro.spike.SpikeFactory.instance(); %Spike Factory (Kaya code)
            sa = sf.getPhyOutputFolder(obj.phypath); %Spike Array (Kaya code)
            sat = SpikeArrayTrack(sa,obj.positionData); %SpikeArrayTrack (Kaya)

            timeWindow = [obj.positionData.time.getStartTimeAbs obj.positionData.time.getEndTimeAbs];
            satWindow = sat.getTimeWindow(timeWindow);

            %[artifactfreeSAT,artifactsaveIdx] = artifactThresholdCircularMaze(obj,sat,timeWindow);
            %[thetathreshSAT,thetasaveIdx] = thetaThresholdCircularMaze(obj,artifactfreeSAT,timeWindow,1);
            
            % unionIdx = union(thetasaveIdx,artifactsaveIdx);
            % 
            % SpikeTimes = satWindow.SpikeTableInSamples.SpikeTimes(unionIdx);
            % SpikeCluster = satWindow.SpikeTableInSamples.SpikeCluster(unionIdx);
            % X = satWindow.SpikeTableInSamples.X(unionIdx);
            % Y = satWindow.SpikeTableInSamples.Y(unionIdx);
            % Z = satWindow.SpikeTableInSamples.Z(unionIdx);
            % 
            % newSat = sat;
            % newSat.SpikeTableInSamples = table(SpikeTimes,SpikeCluster,X,Y,Z);
            
            [sut,acgTask,frm,pfm] = obj.getPlaceFields(sat);

                        classTypes = unique(sa.ClusterInfo.class);
            nclassTypes = length(classTypes);
            nUnits = height(sa.ClusterInfo);


            %%%%%%%%%%%%%%%%%%%%%figure

            for c = 1:nclassTypes
                class = classTypes(c);
                for unit = 1:nUnits %1:length(pfm{c})

                    grid_height = 1; grid_width = 3;
                    h1 = tiledlayout(grid_height,grid_width);
                    title(h1,['Class ' num2str(class) ', unit ' num2str(unit)])
                    position = 1; h = 1; w = 1; ax1 = nexttile(position,[h,w]); %3D plot
                    position = 2; h = 1; w = 1; ax2 = nexttile(position,[h,w]); %place field pair 1
                    position = 3; h = 1; w = 1; ax3 = nexttile(position,[h,w]); %place field pair 1

                    % AX1: plots position data vertically over time, and then the locations of unit firing
                    axes(ax1)
                    try
                        sut{c}{unit}.plotOnTrack3D;
                        hold on
                        sut{c}{unit}.plotPlaneonTrack3D(obj.Events.ruleSwitch)
                        scatter(-50,obj.Event.ruleSwitch,0)
                    catch
                        disp(['Could not plot 3d track for ' num2str(c) ' ' num2str(unit)]);
                    end

                    % AX2: Plots occupancy and firing rate map
                    axes(ax2)
                    try; hold on
                        title('2D position')
                        plot(sut{c}{unit}.TimesInSamples.X,sut{c}{unit}.TimesInSamples.Z,'.')
                        hold off
                    catch
                        disp(['Could not plot 2d position for ' num2str(c) ' ' num2str(unit)]);
                    end

                    % AX2: Plots occupancy and firing rate map
                    axes(ax3)
                    try; hold on
                        title('Firing Rate Map')
                        frm{c}{unit}.plotSmooth;
                        hold off
                    catch;
                        disp(['Could not plot frm for ' num2str(c) ' ' num2str(unit)]);
                    end

                    pause
                    clf
                end
            end

                h = gcf;
                 h.Renderer = 'painters';
                fignamePdf = [obj.basepath '/Figures/3Dplace/position3d_class' num2str(class) '_unit' num2str(unit) '_id' num2str(sa.ClusterInfo.id(unit)) '.pdf'];
                print(fignamePdf,'-dpdf','-r300','-bestfit'); %did smoosh the raster a bit tho
                 fignameFig = [obj.basepath '/Figures/3Dplace/position3d_class' num2str(class) '_unit' num2str(unit) '_id' num2str(sa.ClusterInfo.id(unit)) '.fig'];
                saveas(h,fignameFig)

            %%%%%%%%%%%%%%%%%%%% restart here

            classTypes = unique(sa.ClusterInfo.class);
            nclassTypes = length(classTypes);
            nUnits = height(sa.ClusterInfo);

            %get acg for entire session before cutting down to just task ep
            [~,acg_wholesession] = obj.getEntireSession(sa);
            saTrack = sa.getTimeInterval([obj.Events.taskStart, obj.Events.taskEnd]); 
            [sut_Track,acg_Track] = obj.getEntireSession(saTrack);

            
            %get place fields for individual units
            saTrack_pair1 = sa.getTimeInterval([obj.Events.taskStart, obj.Events.ruleSwitch]); %just for during task for now
            saTrack_pair2 = sa.getTimeInterval([obj.Events.ruleSwitch, obj.Events.taskEnd]); %just for during task for now
           
            [sut1,acgTask1,frm1,pfm1] = obj.getPlaceFields(saTrack_pair1);
            [sut2,acgTask2,frm2,pfm2] = obj.getPlaceFields(saTrack_pair2);

            % plot
            error = [];
            figure
            for c = 1:nclassTypes
                class = classTypes(c);
                for unit = 1:nUnits %1:length(pfm{c})

                    grid_height = 3; grid_width = 5;
                    h1 = tiledlayout(grid_height,grid_width);
                    title(h1,['Class ' num2str(class) ', unit ' num2str(unit)])
                    position = 1; h = 3; w = 1; ax1 = nexttile(position,[h,w]); %3D plot
                    position = 2; h = 2; w = 2; ax2 = nexttile(position,[h,w]); %place field pair 1
                    position = 4; h = 2; w = 2; ax4 = nexttile(position,[h,w]); %place field pair 2
                    position =  12; h = 1; w = 1; ax12 = nexttile(position,[h,w]); % acg whole session
                    position =  13; h = 1; w = 1; ax13 = nexttile(position,[h,w]);  %acg just task
                    position =  14; h = 1; w = 1; ax14 = nexttile(position,[h,w]); % waveform entire session

                    % AX1: plots position data vertically over time, and then the locations of unit firing
                    axes(ax1)
                    %hold on
                    %title('Position Firing over Time')
                    try
                    sut{c}{unit}.plotOnTrack3D;
                    hold on
                    sut{c}{unit}.plotPlaneonTrack3D(obj.Events.ruleSwitch)
                    s%catter(-50,obj.Event.ruleSwitch,0)
                    catch; end
                    %hold off

                    % AX2: Plots occupancy and firing rate map
                    axes(ax2)
                    try; hold on
                    title('Firing Rate Map Pair 1')
                    frm{c}{unit}.plotSmooth;
                    hold off
                    catch; end

                    % AX4: Plots occupancy and firing rate map
                    axes(ax4)
                    try; hold on
                    title('Firing Rate Map Pair 2')
                    frm2{c}{unit}.plotSmooth;
                    hold off
                    catch; end

                    %AX3: Plots acg (for entire session)
                    axes(ax12)
                    try; hold on
                    title('ACG for entire session')
                    acg_wholesession{c}{unit}.plotSingleHistogram
                    hold off
                    catch; end

                      %AX13: Plots acg (for entire session)
                    axes(ax13)
                    try; hold on
                    title('ACG for task epoch')
                    acg_Track{c}{unit}.plotSingleHistogram
                    hold off
                    catch; end

                    %AX4 Plots waveform for entire session

%                     % AX6: Plots place field map % isn't currently working
%                     axes(ax6)
%                     hold on
%                     %title('Place Field Map')
%                     pfm{c}{unit}.plot;
%                     hold off

                    pause
                    clf
                  %  end
%                     catch
%                     error{c}{unit} = 1;
%                     end
            end
            end

        end



        function plotPlaceCellsKD(obj)
            %%%%% uses PFClassic: %%%%%%%
            %pos is nx2 array, NORMALIZED between 0 and 1.
            %SpkCnt gives the number of spikes in each epoch.
            %Smooth is width of gaussian smoother (in 0 ... 1 units)
            % nGrid gives grid spacing (should be larger than 1/smooth)
            %top rate is for the maximum firing rate on the color map, if displayed
            %output is a smoothed occupancy map.

            %%%%%% default vals %%%%%%%% 
            timeWindow = [obj.positionData.time.getStartTimeAbs obj.positionData.time.getEndTimeAbs];          
            posSR = obj.positionData.time.getSampleRate; %position sample rate
            spikeSR = obj.TimeIntervalCombined.getSampleRate;
            lfpSR = obj.LFP.TimeIntervalCombined.getSampleRate;
            smoothVal = 0.02; %out of 1
            nGrid = [100 100]; %one val for each dimension + for time
            timebinsize = .01;% tbins in units of seconds, e.g. 10 ms.

            %%%%% time %%%%%%%%%%%%%%%%
            positionTimeZT = obj.positionData.time.getTimePointsZT;

            % TMP %%%%%%%%%%%%%%%%%%
            %firsthalfindex = 1:floor(length(positionTimeZT)/2);
            %secondhalfindex = (length(firsthalfindex)+1):length(positionTimeZT);
            %positionTimeZT = positionTimeZT(secondhalfindex);
            %%%%%%%%%%%%%%%%%%%%%%%%

            %t = obj.positionData.time.getTimePointsInSamples/posSR; %in seconds
            tbegin = positionTimeZT(1); tend = positionTimeZT(end); 
            nBins = round(seconds(tend-tbegin)/timebinsize);
            %%%%%%%%binned position %%%%%%%%%%
            % Normalize position and get binned position
            x = obj.positionData.data.X;
            z = obj.positionData.data.Z;

            % TMP %%%%%%%%%%%%%%
            %x = x(secondhalfindex);
            %z = z(secondhalfindex);
            %%%%%%%%%%%%%%%%%%%%%%%%%

             x = x - min(x); %make all values positive, aligned to zero
             z = z - min(z);
            % 
             normX = x/max(x); %normalize to make all vals between 0 and 1.
             normZ = z/max(z);

            xzt(:,1) = normX; 
            xzt(:,2) = normZ;
            xzt(:,3) = seconds(positionTimeZT); 
    
            positionBinned = external.Placefields.binpos(xzt,nBins);
            histedges = positionBinned(:,3)-timebinsize;
            histedges(end+1) = histedges(end)+timebinsize;
                     
            % Set limits defined by the shape of the track
           % limits.x = [0.3 0.7];
           % limits.z = [0.3 0.7];

            % get eeg states
           % eegStates = getStates(obj,timeWindow);
            basepath = [obj.basepath '/'];
           % rejectCh = [97:127];
           % thetaCh = [64:74];
            % SleepState = SleepScoreMaster_RW(obj.basepath,'Notch60Hz',1,'overwrite',true);%,'rejectChannels',rejectCh,'thetaChannels',thetaCh,'overwrite',
            % Get theta power at each point in 

            [thetaPow,thetaTfm] = getPowerTimeWindow(obj,timeWindow,[150 200]);
            lfpTimeInterval = obj.LFP.TimeIntervalCombined.getTimeIntervalForTimes(timeWindow);
            lfpTime = seconds(lfpTimeInterval.getTimePointsZT);   
            
            thetaBinned = interp1(lfpTime,thetaPow.Values,positionBinned(:,3)); %find theta power at the times where lfpTime = binnedPos time.

            thetaThresh = 1;
            zTheta = (thetaBinned- nanmean(thetaBinned))/nanstd(thetaBinned);
            %figure
             %thetaHist = histogram(zTheta,100);
            belowthreshTheta = find(zTheta < 0); %ARBITRARY IN MIDDLE RIGHT NOW, no bimodality

            %%%%%% spiking info %%%%%%%%%%
            sf = neuro.spike.SpikeFactory.instance(); %Spike Factory (Kaya code)
            sa = sf.getPhyOutputFolder(obj.phypath); %Spike Array (Kaya code)
            su = sa.getSpikeUnits; % pull spiketimes per unit.
            nUnits = length(su);

            classTypes = unique(sa.ClusterInfo.class);
            nclassTypes = length(classTypes);            

            %sat = neuro.spike.SpikeArrayTrack(sa,obj.positionData); %SpikeArrayTrack (Kaya)
            %sut = obj.getPlaceFields(sat); %get just the spike info for the track period

            %%%%%% for each unit, plot place field %%%%%%%
            figure
            grid_height = 2; grid_width = 3;
            for c = 1 %:nclassTypes
              %  class = classTypes(c);
                for unit = 1:nUnits %1:length(pfm{c})

                     %%%%%% get spikeCount %%%%%%%%
                    % spikeCount is binned with same bins as binnedPos above
                    
                    spksZT = seconds(su(unit).getTimesZT)';  %spikes in ZT time.                

                    spksBinned = histcounts(spksZT,histedges)';
                    spksBinned(end)=0; spksBinned(1)=0; % find number of spikes per unit time.
                  
                    spkslogical = logical(spksBinned);
                    xprevelocityfilter = positionBinned(spkslogical,1)*100;
                    zprevelocityfilter = positionBinned(spkslogical,2)*100;

                    %%%%%%%%% % velocityfilter = external.velocityfilter(binnedPos)
                    velocityThresh = 10; %in 10 cm/hour
                    [~,~,spksBinned] = external.velocityfilterUKRW_2D(positionBinned,velocityThresh,spksBinned); %in cm/s
                   % spkCount2(belowthreshTheta) = 0;
                  %  [xztRange, nspkRange] = position.SetCircularRange(binnedPos,spkCount2,limits);

                    spkslogical = logical(spksBinned);
                    xpostvelocityfilter = positionBinned(spkslogical,1)*100;
                    zpostvelocityfilter = positionBinned(spkslogical,2)*100;

                    h1 = tiledlayout(grid_height,grid_width);
                    title(h1,['Unit ' num2str(unit)]);
                    position = 1; h = 1; w = 1; ax1 = nexttile(position,[h,w]); %3D plot
                    position = 2; h = 1; w = 1; ax2 = nexttile(position,[h,w]); %place field pair 1
                    %position = 3; h = 1; w = 1; ax3 = nexttile(position,[h,w]); %place field pair 1
                    position = 4; h = 1; w = 1; ax4 = nexttile(position,[h,w]); %place field pair 1
                    position = 5; h = 1; w = 1; ax5 = nexttile(position,[h,w]); %place field pair 1
                    position = 6; h = 1; w = 1; ax6 = nexttile(position,[h,w]); %place field pair 1

                    % AX1: plots position data vertically over time, and then the locations of unit firing
                    axes(ax1)
                    hold on
                    title('Spk locs. (all velocities)')
                    plot(x,z,'y')
                    plot(xprevelocityfilter,zprevelocityfilter,'.r')
                    %plot(sut{1}{unit}.TimesInSamples.X,sut{1}{unit}.TimesInSamples.Z,'.r');
                    ax1.DataAspectRatio = [1 1 1];
                    hold off

                    % AX2: plots position data vertically over time, and then the locations of unit firing
                    axes(ax2)
                    hold on
                    title('Spk locs. (> 10 cm velocities)')
                    plot(x,z,'y')
                    plot(xpostvelocityfilter,zpostvelocityfilter,'.r')
                    ax2.DataAspectRatio = [1 1 1];
                    %plot(sut{1}{unit}.TimesInSamples.X,sut{1}{unit}.TimesInSamples.Z,'.r');
                    hold off


                    % AX3: Plots occupancy and firing rate map
                    axes(ax6)
                    hold on

                    %Get and plot place field %%%%%%%
                    [~,~,nSpikes,sTimeSpent,snSpikes] = ...
                        external.Placefields.PFClassic(positionBinned(:,1:2),spksBinned,smoothVal,nGrid,timebinsize);
                       
                     title(['Place map, nspikes = ' num2str(sum(sum(nSpikes)))]) 

                    %   linkaxes([ax1 ax2],'xy')

                    % AX4: 
                    axes(ax4)
                    hold on
                    title('smoothed Occ. map')
                    imagesc(sTimeSpent')
                    ylim([0 length(sTimeSpent)])
                    xlim([0 length(sTimeSpent)])
                    ax4.DataAspectRatio = [1 1 1];
                    hold off

                    % AX5: 
                    axes(ax5)
                    hold on
                    title('smoothed Spikemap')
                    imagesc(snSpikes')
                    ylim([0 length(snSpikes)])
                    xlim([0 length(snSpikes)])
                    ax5.DataAspectRatio = [1 1 1];
                    hold off

                    pause
                    clf

                end
            end

        end

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

        function [sut,acgTask,frm,pfm] = getPlaceFields(obj,sa)
            spikeUnits = sa.getSpikeUnits();
            classTypes = 1; %just for UK data
            %classTypes = unique(sa.ClusterInfo.class);

            for c = 1%:length(classTypes)
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






