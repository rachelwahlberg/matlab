classdef (Abstract) FigurePipeline < handle

    %%FigurePipeline %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Rachel Wahlberg Aug 2023 %
    % This pipeline is largely based off of Utku Kaya code %

    %      basepath 
    %         basename 
    %         behaviorpath
    %         phyoutputpath =
    %         TimeIntervalCombinedpath = ...
    %             '/data/ExperimentsRSW/CircularMaze/20230614/merged_20230614/merged_20230614'
    %         fullfile(basepath, [basename '_1250Hz.TimeIntervalCombined.csv']);
    %         phyfolder = '/home/wahlberg/merged_20230614crs.GUI/';
    %         optifolder = 
    % end

    % the < handle fixed some property inheritance issues

    %% SETUP

    properties %constant means subclasses can't change these properties
        basepath = '/data/ExperimentsRSW/CircularMaze/Harry/CA1/20230614/merged_20230614';
        basename = 'merged_20230614';
        behaviorpath = '/data/ExperimentsRSW/CircularMaze/Harry/CA1/20230614/merged_20230614/behaviorfiles/';
        phypath = '/home/wahlberg/merged_20230614crs.GUI/';
        optipath = '/data/ExperimentsRSW/CircularMaze/Harry/CA1/20230614/Optitrack';
        SpikeTableInSamples
        ClusterInfo
        Info
        positionData
        TimeIntervalCombined
        Events
    end % does NOT include LFP for more flexibility on channel selection

    methods(Abstract)
        testSubplots(obj)
        plotFigure(obj)
        saveFigure(obj)
    end

    %% FUNCTIONS

    methods

        function initializeProperties(obj)

            cd(obj.basepath)
            %%%%%%%%%%%%%%%%% spikes %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            sf = neuro.spike.SpikeFactory.instance(); %Spike Factory (Kaya code)
            sa = sf.getPhyOutputFolder(obj.phypath); %Spike Array (Kaya code)

            obj.SpikeTableInSamples = sa.SpikeTableInSamples;
            obj.ClusterInfo = sa.ClusterInfo;
            obj.Info = obj.Info;
            
            %get absoulute times for epoch start/end and rule switch
            obj.Events = readtable([obj.basepath '/' obj.basename '_events.xlsx']); %xlsx allows format to stay, whereas csv wouldn't

           % p = obj.loadPositionData; % a position.PositionData object.
            obj.positionData = obj.loadPositionData;

            %%%%%%%%%%%%%% timestamps %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            theFile=dir(fullfile(obj.phypath,['*TimeIntervalCombined*' '.csv'])); %TIC file gives timing info for each pre-merge file
            if isempty(theFile)
                logger.warning(strcat('TimeIntervalCombined is not loaded. \n\tLocation:\t',foldername,'\n'));
            end

            ticd=time.TimeIntervalCombined(fullfile(theFile.folder, theFile.name)); % from the time folder
            ticd = ticd.setZeitgeberTime(hours(12)); % for noon start, dark to light transition
            obj.TimeIntervalCombined = ticd; %pulling the 1250Hz version!

        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%% LOAD POSITION DATA %%%%%%%%%%%%%%%%%%

        function positionData = loadPositionData(obj)
            fnames = dir(fullfile(obj.optipath,sprintf('*Take*')));
            oploader = position.optiTrack.OptiLoader.instance();
            %filenames =  {fnames(:).name};
            for f = 1:length(fnames)
                filename= fullfile(obj.optipath,fnames(f).name);
                op = oploader.loadFile(filename);
                pos0 = op.Positions;
                savefolder = [];
                pd_unmerged{f}=position.fillRigidBodyNans(pos0,savefolder);

                if ~exist("positionData")
                    positionData = pd_unmerged{f};
                else
                    positionData.data = [positionData.data; pd_unmerged{f}.data(:,{'X','Y','Z'})];
                    positionData.time = positionData.time + pd_unmerged{f}.time;
                end
            end
      end

%         function interpolatedData = interpPositionData(obj) %OLD (now use fillRigidBodyNans)
% 
%             interpolatedData = obj.positionData;
%             rawx = obj.positionData.data.X; 
%             rawy = obj.positionData.data.Y; 
%             rawz = obj.positionData.data.Z; 
% 
%             %rewrite
%             interpolatedData.data.X = interp1(t,rawx,spks,'nearest');   % Align the position of the animal when the cell spikes to an x-grid position
%             yint = interp1(t,yb,spks,'nearest');   % Align the position of the animal when the cell spikes to an y-grid position
%         end

        function [] = getWavelet(obj,timeframe,freq)
            % timeframe in [start end] format, in seconds
            % frequency in hz, [low high]
            %%%%%%%%%%% Buzcode wavelet funtion (Morlet) %%%%%%%%%%%%%%%%%%%%%%%%%%
            %check where the ripple power is

            whitened = WhitenSignal(obj.LFP.Values);
            sr = obj.TimeIntervalCombined.getSampleRate;

            lfpSegment = whitened(floor(timeframe(1)*sr):ceil(timeframe(2)*sr));

            wavespec = bz_WaveSpec(lfpSegment,'samplingRate',sr,'frange',[freq(1) freq(2)],...
                'nfreqs',300,'showprogress',true,'ncyc',7);
            % trade off number of frequencies and range of frequencies in order to not overwhelm memory

            figure
            hold all
            title(['Ch ' num2str(obj.LFP.ChannelName) ', ' ...
                num2str(timeframe(1)) '-' num2str(timeframe(2)) ' sec' ]);
            imagesc(wavespec.timestamps,wavespec.freqs,(abs(wavespec.data))');
            ax = gca;
            ax.YDir = 'normal'; ax.YScale="log"; ax.YLim = [0.5 300]; ax.XLim = [min(wavespec.timestamps) max(wavespec.timestamps)];
            colorbar; ax.CLim=[0 300];
        end

%         %%%%%%%%%%%%%%%%%%%%%%%%%% CREATE METHODS FIGURE %%%%%%%%%%%%%%%%%%
%         function [] = callMethods()
%             filepathtosave = [basepath '/Figures/Methods/methods_' date '.fig'];
%             createMethodsFigure(filepathtosave);
%         end
% 
%         %%%%%%%%%%%%%%%%%%% CREATE ONLINE PLACE CELLS FIGURE %%%%%%%%%%%%%%
%         function [] = callOnlinePlaceCells
%             % a) ex place cell at beginning of task, end of pair 1, end of pair 2
%             % b) heat map place cell as a function of time
%             % c) fidelity?
%             % d) jumpiness?
%             % e) etc...
%         end

        %%%%%%%%%%%%%%%%%%%%% CREATE ONLINE THETA FIGURE %%%%%%%%%%%%%%%%%%
% 
%         function [] = callOnlineTheta
%             % a) th power as func of time
%             % b) speed as func of time
%             % c) example locked raster over time
%             % d) performance as func of time
%             % e) summary stats...
%             % f) more summary stats ...
%             % g) look ahead distance
%             % h) other redish comparisons?
%         end

        %%%%%%%%%%%%%%%%%%%%%% CALL SWR INTRO FIGURE %%%%%%%%%%%%%%%%%%%%%%

%         function [] = callSWRintro
%             % a) SWR power vs #replays
%             % b) example replay
%             % c) ripple power vs performance
%             % d) etc...
%         end

        %%%%%%%%%%%%%%%%%% CALL REPLAY OVER TIME FIGURE %%%%%%%%%%%%%%%%%%%
% 
%         function [] = callReplayOverTime
%             % a) examples over time
%             % b) ask pho for inspo (with his long vs short track stuff)
%         end
% 
%         %%%%%%%%%%%%%%%%%%%% CALL COLOR ANALYSIS FIGURE %%%%%%%%%%%%%%%%%%%
% 
%         function [] = callColorAnalyses
%             % a) task rule (drawing)
%             % b) example place cell at boundary
%             % c) violin plots of #PC at boundary in cross versus same color pairs
%             % d) black performance vs white vs cross
%             % e) time to performance for b vs w vs c
%             % f) etc...
%         end
% 
%         %%%%%%%%%%%%%%%%%%%% CALL DENTATE GYRUS FIGURE %%%%%%%%%%%%%%%%%%%%
% 
%         function [] = callDentateGyrusAnalyses
%             % a) example cells
%             % b) violin plots summary stats
%             % c) PSD
%             % d) raster
%             % e) etc... (go through hpc analyses but with dg cells)
% 
%         end

    end

end