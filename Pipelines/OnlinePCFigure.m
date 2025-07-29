classdef OnlinePCFigure < FigurePipeline
% Collection of figures and analyses for online place cell analysis (no
% replay, etc)

%%%%% METHODS %%%%%%%%%%%%%%%%%%%%%%%%%%
% initializeProperties (inherited from FigurePipeline)
% loadPositionData (inherited from FigurePipeline)
% testSubplots
% plotFigure
% saveFigure




    properties
        LFP
    end

    methods

        function obj = OnlinePCFigure()
            % to call, onPC = OnlinePCFigure()

            obj.initializeProperties(); % also cd's to basepath

            %%%%%%%%%%%%%%%%% lfp %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            ctdh = neuro.basic.ChannelTimeDataHard(obj.basepath);
            lfpchannel = 30;
            obj.LFP = ctdh.getChannel(lfpchannel);
            ticd = obj.LFP.TimeIntervalCombined; 
            ticd = ticd.setZeitgeberTime(hours(12)); % for noon start, dark to light transition
            obj.LFP.TimeIntervalCombined = ticd; 

            % get online times - no rest periods

        end

        function [] = plotPhasePrecession(obj)
        
            %1) get all spiketimes for a unit
            %2) bandpass lfp in theta band
            %3) get phase of theta for each spiketime
            %4) sort by theta phase
            %5) plot theta phase against time for each unit

            % each place cell in a cell
            %idx = 'class';
            import neuro.spike.*
            import neuro.placeField.*

            % pd_ontrack = obj.positionData.getTimeWindow(timewindow)
            sf = neuro.spike.SpikeFactory.instance(); %Spike Factory (Kaya code)
            sa = sf.getPhyOutputFolder(obj.phypath); %Spike Array (Kaya code)

            spikeUnits = sa.getSpikeUnits;  

            figure

            for iU = 1:length(spikeUnits)
                class = num2str(spikeUnits(iU).Info.class);
                unit = num2str(spikeUnits(iU).Info.id);

                if strcmp(class,"1") ~= 1 && strcmp(class,"2") ~= 1
                    continue
                end

            ztTimes = spikeUnits(iU).getTimesZT;
            spikeUnitsClockwise = spikeUnits(iU);
                  spikeUnitsClockwise.TimesInSamples = [];
             
                cl_counter = 1;
                for s = 1:length(ztTimes)
                    cl = abs(ztTimes(s)-obj.Events.Clockwise);
                    if min(cl)<duration('00:00:00:001')
                        spikeUnitsClockwise.TimesInSamples.SpikeTimes(cl_counter) = spikeUnits(iU).TimesInSamples.SpikeTimes(s);
                        cl_counter = cl_counter + 1;
                    end
                end
        
                hold on
                h1 = tiledlayout(1,1);
                
            %  title(h1,['Class ' class ', unit ' unit])
               position = 1; h = 1; w = 1; ax1 = nexttile(position,[h,w]); %3D plot
              
                %clockwise object
                try
                    
                    sulCl = SpikeUnitLFP(spikeUnitsClockwise(iU), obj.LFP);
                             
                    spikePhases = sulCl.getPhases; %a spikePhases class
                    
                    absTimes = sulCl.getAbsoluteSpikeTimes;
                    stSampleRate = sulCl.Time.getSampleRate;
                    posDataClockwise = obj.positionData.getPositionForTimes(absTimes');%,stSampleRate);
                    pp = neuro.phase.PhasePrecession(spikePhases,posDataClockwise,ztTimes');

                    axes(ax1)
                    pp.plotHist;

                catch
                    continue
                end

                pause

                h = gcf;
                h.Renderer = 'painters';
           %     figname = [obj.basepath '/Figures/ThetaPhase/thetaphase_badprecession_class' class '_unit' num2str(iU) '_id' unit '.pdf'];
               figname = [mf.basepath '/Figures/behavior_4thyeartalk.pdf'];
                print(figname,'-dpdf','-r300','-bestfit'); %did smoosh the raster a bit tho

                clf
                hold off
            end
        end
 



        function [] = plotPhasePrecessionDirections(obj)
        
            %1) get all spiketimes for a unit
            %2) bandpass lfp in theta band
            %3) get phase of theta for each spiketime
            %4) sort by theta phase
            %5) plot theta phase against time for each unit

            % each place cell in a cell
            %idx = 'class';
            import neuro.spike.*
            import neuro.placeField.*

            % pd_ontrack = obj.positionData.getTimeWindow(timewindow)
            sf = neuro.spike.SpikeFactory.instance(); %Spike Factory (Kaya code)
            sa = sf.getPhyOutputFolder(obj.phypath); %Spike Array (Kaya code)


            spikeUnits = sa.getSpikeUnits;  

            figure

            for iU = 1:length(spikeUnits)
                class = num2str(spikeUnits(iU).Info.class);
                unit = num2str(spikeUnits(iU).Info.id);

                if strcmp(class,"1") ~= 1 && strcmp(class,"2") ~= 1
                    continue
                end

                ztTimes = spikeUnits(iU).getTimesZT;
                spikeUnitsClockwise = spikeUnits(iU);
                spikeUnitsCounter = spikeUnits(iU);
                spikeUnitsClockwise.TimesInSamples = [];
                spikeUnitsCounter.TimesInSamples = [];

                cl_counter = 1;
                co_counter = 1;
                for s = 1:length(ztTimes)
                    cl = abs(ztTimes(s)-obj.Events.Clockwise);
                    if min(cl)<duration('00:00:00:001')
                        spikeUnitsClockwise.TimesInSamples.SpikeTimes(cl_counter) = spikeUnits(iU).TimesInSamples.SpikeTimes(s);
                        cl_counter = cl_counter + 1;
                    else 
                        co = abs(ztTimes(s)-obj.Events.Counter);
                        if min(co)<duration('00:00:00:001')
                        spikeUnitsCounter.TimesInSamples.SpikeTimes(co_counter) = spikeUnits(iU).TimesInSamples.SpikeTimes(s);
                        co_counter = co_counter + 1;
                        end
                    end
                end

                grid_height = 2; grid_width = 2;
                h1 = tiledlayout(grid_height,grid_width);
                title(h1,['Class ' num2str(class) ', unit ' num2str(unit)])
                position = 1; h = 1; w = 1; ax1 = nexttile(position,[h,w]); %3D plot
                position = 2; h = 1; w = 1; ax2 = nexttile(position,[h,w]); %place field pair 1
                position = 3; h = 1; w = 1; ax3 = nexttile(position,[h,w]); %place field pair 2
                position = 4; h = 1; w = 1; ax4 = nexttile(position,[h,w]); % acg whole session

                %clockwise object
                try
                sulCl = SpikeUnitLFP(spikeUnitsClockwise, obj.LFP);
                spikePhasesCl = sulCl.getPhases; %a spikePhases class
                ztTimesCl = sulCl.getTimesZT;
                absTimesCl = sulCl.getAbsoluteSpikeTimes;
                stSampleRate = sulCl.Time.getSampleRate;
                posDataClockwise = obj.positionData.getPositionForTimes(absTimesCl');%,stSampleRate);
                ppClockwise = neuro.phase.PhasePrecession(spikePhasesCl,posDataClockwise,ztTimesCl');
                
                axes(ax1) %clockwise trials
                ppClockwise.plotPrecession; %theta phase against position
                
                axes(ax3) %phase distribution circular histogram
                ppClockwise.plotHist;
                
                catch
                end

                try
                % counterclockwise object
                sulCo = SpikeUnitLFP(spikeUnitsCounter, obj.LFP);
                spikePhasesCo = sulCo.getPhases; %a spikePhases class
                ztTimes = sulCo.getTimesZT;
                posDataCounter = obj.positionData.getPositions_largeArrays(ztTimes',stSampleRate);
                ppCounter = neuro.phase.PhasePrecession(spikePhasesCo,posDataCounter,ztTimes');
                
                axes(ax2)
                ppCounter.plotPrecession; %

                axes(ax4)
                ppCounter.plotHist;

                catch
                end
                
                pause
                clf
            end
        end

        function [clockwiseStats,counterStats] = getStatsDirection(obj)
        
            %1) get all spiketimes for a unit
            %2) bandpass lfp in theta band
            %3) get phase of theta for each spiketime
            %4) sort by theta phase
            %5) plot theta phase against time for each unit

            % each place cell in a cell
            %idx = 'class';
            import neuro.spike.*
            import neuro.placeField.*

            % pd_ontrack = obj.positionData.getTimeWindow(timewindow)
            sf = neuro.spike.SpikeFactory.instance(); %Spike Factory (Kaya code)
            sa = sf.getPhyOutputFolder(obj.phypath); %Spike Array (Kaya code)
            spikeUnits = sa.getSpikeUnits;  

            tmp = nan(length(spikeUnits),14);

            clockwiseStats = array2table(tmp,...
                'VariableNames',{'Mean' 'Median' 'Var' 'Std' 'Std0' 'Skewness' ...
                'Skewness0' 'Kurtosis' 'Kurtosis0' 'RayleighP' 'RayleighZ' ...
                'OmnibusP' 'OmnibusZ' 'SymmetryP'});
            counterStats = array2table(tmp,...
                'VariableNames',{'Mean' 'Median' 'Var' 'Std' 'Std0' 'Skewness' ...
                'Skewness0' 'Kurtosis' 'Kurtosis0' 'RayleighP' 'RayleighZ' ...
                'OmnibusP' 'OmnibusZ' 'SymmetryP'});

            for iU = 1:length(spikeUnits)
                class = num2str(spikeUnits(iU).Info.class);
                unit = num2str(spikeUnits(iU).Info.id);

                if strcmp(class,"1") ~= 1 && strcmp(class,"2") ~= 1
                    continue
                end

                ztTimes = spikeUnits(iU).getTimesZT;
                spikeUnitsClockwise = spikeUnits(iU);
                spikeUnitsCounter = spikeUnits(iU);
                spikeUnitsClockwise.TimesInSamples = [];
                spikeUnitsCounter.TimesInSamples = [];

                cl_counter = 1;
                co_counter = 1;
                for s = 1:length(ztTimes)
                    cl = abs(ztTimes(s)-obj.Events.Clockwise);
                    if min(cl)<duration('00:00:00:001')
                        spikeUnitsClockwise.TimesInSamples.SpikeTimes(cl_counter) = spikeUnits(iU).TimesInSamples.SpikeTimes(s);
                        cl_counter = cl_counter + 1;
                    else 
                        co = abs(ztTimes(s)-obj.Events.Counter);
                        if min(co)<duration('00:00:00:001')
                        spikeUnitsCounter.TimesInSamples.SpikeTimes(co_counter) = spikeUnits(iU).TimesInSamples.SpikeTimes(s);
                        co_counter = co_counter + 1;
                        end
                    end
                end

                try
                sulCl = SpikeUnitLFP(spikeUnitsClockwise, obj.LFP);
                spikePhasesCl = sulCl.getPhases; %a spikePhases class             

                tmpClockStats = spikePhasesCl.getStats;
                clockwiseStats.Mean(iU) = tmpClockStats.mean;
                clockwiseStats.Median(iU) = tmpClockStats.median;
                clockwiseStats.Var(iU) = tmpClockStats.var;
                clockwiseStats.Std(iU) = tmpClockStats.std;
                clockwiseStats.Std0(iU) = tmpClockStats.std0;
                clockwiseStats.Skewness(iU) = tmpClockStats.skewness;
                clockwiseStats.Skewness0(iU) = tmpClockStats.skewness0;
                clockwiseStats.Kurtosis(iU) = tmpClockStats.kurtosis;
                clockwiseStats.Kurtosis0(iU) = tmpClockStats.kurtosis0;
                [clockwiseStats.RayleighP(iU),clockwiseStats.RayleighZ(iU)] = spikePhasesCl.getTestRayleigh;
                [clockwiseStats.OmnibusP(iU),clockwiseStats.OmnibusZ(iU)] = spikePhasesCl.getTestOmnibus;
                clockwiseStats.SymmetryP(iU) = spikePhasesCl.getTestSymmetryAroundMedian;
                catch
                end

                try
                sulCo = SpikeUnitLFP(spikeUnitsCounter, obj.LFP);
                spikePhasesCo = sulCo.getPhases; %a spikePhases class

                tmpCounterStats = spikePhasesCo.getStats;
                counterStats.Mean(iU) = tmpCounterStats.mean;
                counterStats.Median(iU) = tmpCounterStats.median;
                counterStats.Var(iU) = tmpCounterStats.var;
                counterStats.Std(iU) = tmpCounterStats.std;
                counterStats.Std0(iU) = tmpCounterStats.std0;
                counterStats.Skewness(iU) = tmpCounterStats.skewness;
                counterStats.Skewness0(iU) = tmpCounterStats.skewness0;
                counterStats.Kurtosis(iU) = tmpCounterStats.kurtosis;
                counterStats.Kurtosis0(iU) = tmpCounterStats.kurtosis0;
                [counterStats.RayleighP(iU),counterStats.RayleighZ(iU)] = spikePhasesCo.getTestRayleigh;
                [counterStats.OmnibusP(iU),counterStats.OmnibusZ(iU)] = spikePhasesCo.getTestOmnibus;
                counterStats.SymmetryP(iU) = spikePhasesCo.getTestSymmetryAroundMedian;
                catch
                end
                    iU
            end
        end


        function [] = testSubplots(obj)

        end

        function [] = plotFigure(obj)

        end

        function [] = saveFigure(obj)

        end
    end

end