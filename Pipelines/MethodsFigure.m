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

        function [] = getSleepStates(obj)

            %%%%%%%%%%%%%Check EEG states with sliding window %%%%%%%%%%%%%%%%%%%%%%%%%
            % Put "100" as a window size into WhitenSignal, didn't change much
            % figure out what the measurement is for the window

            %%%%%%%%%%%%%Check EEG states for entire behavioral period %%%%%%%%%%%%%%%%
            CheckEegStates_RW(obj.basepath,obj.basename,'lfp',(obj.LFP.Values)','redo_flag',true); % and testpre, and testpost
            %select ch 1 because I'm just sending in 1 channel anyway

            %%%%%%%%%%%%%% Try SleepScoreMaster for theta detection %%%%%%%%%%%%%%%%%%

            %badCh = session.channelTags.Bad.channels;
            %thetaCh = session.channelTags.ThetaChan.channels; % check this

%             badCh = [0 122];
%             thetaCh = 126;
%             SleepState = SleepScoreMaster(basepath,'ignoretime',removetimestamps.sampformat,'rejectchannels',badCh,'ThetaChannels',thetaCh);
         end

        function [] = getPSD(obj)


        end

        function [] = spikeRaster(obj)

            %%first argument is what column to sort into different figures, second arg
            %is the filepath to which to save figures (kills matlab with too many
            %units!)
            sf = neuro.spike.SpikeFactory.instance(); %Spike Factory (Kaya code)
            sa = sf.getPhyOutputFolder(obj.phypath); %Spike Array (Kaya code)

            % is currently creating a figure of its own
            sa.plotRaster('sh','/home/wahlberg/merged_20230614crs.GUI/Figures/')
        end


        function [] = plotPlaceMaps(obj,idx) %Occupancy currently looks cut off!
    
            % each place cell in a cell
            %idx = 'class';
            import neuro.spike.*
            import neuro.placeField.*

            % pd_ontrack = obj.positionData.getTimeWindow(timewindow)
            sf = neuro.spike.SpikeFactory.instance(); %Spike Factory (Kaya code)
            sa = sf.getPhyOutputFolder(obj.phypath); %Spike Array (Kaya code)
            sat = SpikeArrayTrack(sa,obj.positionData); %SpikeArrayTrack (Kaya)
            
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
            for c = 1%:nclassTypes
                class = classTypes(c);
                for unit = 1%:nUnits %1:length(pfm{c})

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
                    sut_Track{c}{unit}.plotOnTrack3D;
                    hold on
                    sut_Track{c}{unit}.plot3DEvent(obj.Events)
                    s%catter(-50,obj.Event.ruleSwitch,0)
                    catch; end
                    %hold off

                    % AX2: Plots occupancy and firing rate map
                    axes(ax2)
                    try; hold on
                    title('Firing Rate Map Pair 1')
                    frm1{c}{unit}.plotSmooth;
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

        function [sut,acg] = getEntireSession(obj,sa)

            % get ACGs and waveforms BEFORE grabbing out just the task ep
            spikeUnits = sa.getSpikeUnits();
            classTypes = unique(sa.ClusterInfo.class);
            for c = 1%:3%length(classTypes)
                class = classTypes(c);
                unitcount = 1;
                for unit = 1%:length(spikeUnits)
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
            classTypes = unique(sa.ClusterInfo.class);

            for c = 1%:3%length(classTypes)
                class = classTypes(c);
                unitcount = 1;
                for unit = 1%:length(spikeUnits)
                    if spikeUnits(unit).Info.class == class
                        try
                            sut{c}{unitcount} = neuro.spike.SpikeUnitTracked(spikeUnits(unit),obj.positionData);
                            acgTask{c}{unitcount} = sut{c}{unitcount}.getACC;
                            frm{c}{unitcount} = sut{class}{unitcount}.getFireRateMap; %underlay is occupancy map
                            pfm{c}{unitcount} = frm{c}{unitcount}.getPlaceFieldMap;
                        catch
                            sut{c}{unitcount} = [];
                            acgTask{c}{unitcount} = [];
                            pfm{c}{unitcount} = [];
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





