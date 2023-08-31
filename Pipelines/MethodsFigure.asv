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


        function placeCells = getPlaceCells(obj)
            % each place cell in a cell
            speedthreshold = 5;

            sf = neuro.spike.SpikeFactory.instance(); %Spike Factory (Kaya code)
            sa = sf.getPhyOutputFolder(obj.phypath); %Spike Array (Kaya code)
            spikeUnits = sa.getSpikeUnits();
                
            spikeUnits.plotPlaceFieldBoth(obj.positionData,speedthreshold)
            
            







            neuro.placecells.tuningcurves2D(behaviortimestamps,phyoutputpath,'trials',trials,'labelstouse',1:2,'showplots',true); %false);

            [FRmaps,xbins,ybins,empties] = place_cells03_RW(behaviortimestamps,optitrack,phyoutputpath,'labelstouse',1:2,'placecellplots',true); %false);

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





