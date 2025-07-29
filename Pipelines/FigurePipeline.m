classdef (Abstract) FigurePipeline < handle

    %%FigurePipeline %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Rachel Wahlberg Aug 2023 %
    % This pipeline is largely built from Utku Kaya code %
    % the < handle fixed some property inheritance issues
    
    %%%%%%%%%%% SUBCLASSES: %%%%%%%%%%%%%%%%%
    % MethodsFigure
    % OnlinePlaceCellsFigure

    %%%%%%%%%%% METHODS %%%%%%%%%%%%%%%%%%%%%
    %initializeProperties
    %loadPositionData
    %testSubplots (abstract)
    %plotFigure(abstract)
    %saveFigure(abstract)


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

    

    %% SETUP

    properties %constant means subclasses can't change these properties
        basepath = '/data/ExperimentsRSW/CircularMaze/Harry/CA1/20230614/merged_20230614_2';
        basename = 'merged_20230614_2';
        behaviorpath = '/data/ExperimentsRSW/CircularMaze/Harry/CA1/20230614/merged_20230614_2/behaviorfiles/';
        phypath = '/data/ExperimentsRSW/CircularMaze/Harry/CA1/20230614/merged_20230614_2/merged_20230614_2/merged_20230614_2crs.GUI/';
        optipath = '/data/ExperimentsRSW/CircularMaze/Harry/CA1/20230614/Optitrack';
        badChannels = [96:127];
        thetaChannel = 70;
        EMGChannels = [1 31 64 95]; % 4 inmost of each working shank
        % basepath = '/data/UK_ExampleData/AG_2019-12-26_SD';
        % basename = 'AG_2019-12-26_SD';
        % phypath = '/data/UK_ExampleData/AG_2019-12-26_SD/Units/s2';
        % optipath = '/data/UK_ExampleData/AG_2019-12-26_SD/_Position';
        % badChannels = [];
        % thetaChannel = 15;
        % EMGChannels = [0 23 24 25];
        circleParams
        SpikeTableInSamples
        ClusterInfo
        Info
        positionData
        TimeIntervalCombined
        Events
        positionAngle
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

            %get ABSOLUTE times for epoch start/end and rule switch
            obj.Events = cell([]);
            obj.Events.taskTimestamps = readtable([obj.basepath '/' obj.basename '_events.xlsx']); %xlsx allows format to stay, whereas csv wouldn't

            %obj.positionData = obj.loadPositionData('/data/UK_ExampleData/AG_2019-12-26_SD/_Position/Take 2019-12-26 12.54.59 PM_TRACK.csv');
            obj.positionData = obj.loadPositionData;
            %%%%%%%%%%%%%% timestamps %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            theFile=dir(fullfile(obj.phypath,['*TimeIntervalCombined*' '.csv'])); %TIC file gives timing info for each pre-merge file
            %theFile=dir(fullfile(obj.phypath,'..',['*TimeIntervalCombined*' '.csv'])); %UK VERSION
            
            if isempty(theFile)
                logger.warning(strcat('TimeIntervalCombined is not loaded. \n\tLocation:\t',foldername,'\n'));
            end

            ticd=time.TimeIntervalCombined(fullfile(theFile.folder, theFile.name)); % from the time folder
            ticd = ticd.setZeitgeberTime(hours(12)); % for noon start, dark to light transition
            obj.TimeIntervalCombined = ticd; %pulling the 1250Hz version!

            obj.circleParams.innerRadius = 35;
            obj.circleParams.outerRadius = 41;
            obj.circleParams.circCenter = [48,51];
            obj.circleParams.trackRadius = 38;
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%% POSITION FUNCTIONS %%%%%%%%%%%%%%%%%%

        function positionData = loadPositionData(obj,usefile)
            if exist("usefile","var")
                fnames = dir(usefile); %have the full pathname
            else
            fnames = dir(fullfile(obj.optipath,sprintf('*Take*')));
            end
            oploader = position.optiTrack.OptiLoader.instance();
            %filenames =  {fnames(:).name};
            for f = 1:length(fnames)
                filename= fullfile(obj.optipath,fnames(f).name);
                op = oploader.loadFile(filename);
                pos0 = op.Positions;
                savefolder = [];
                pd_unmerged{f}=position.fillRigidBodyNans(pos0,savefolder);

               %%% this is for if you want just the rigid body %%%
                % rigidbody_idx = find(strcmp(op.Positions{:,1},'Rigid Body') == 1);
                % pd_unmerged{f}=op.Positions{rigidbody_idx,5}{1};
                % 
                % pd_unmerged{f}.data.X = pd_unmerged{f}.data.X * 23;
                % pd_unmerged{f}.data.Y = pd_unmerged{f}.data.Y * 23;
                % pd_unmerged{f}.data.Z = pd_unmerged{f}.data.Z * 23;
                % idx1=pd_unmerged{f}.data.X<-60|pd_unmerged{f}.data.X>73;
                % idx2=pd_unmerged{f}.data.Y<-20|pd_unmerged{f}.data.Y>34;
                % idx3=pd_unmerged{f}.data.Z<-10|pd_unmerged{f}.data.Z>121;
                % idxall=idx1|idx2|idx3;
                % pd_unmerged{f}.data.X(idxall)=nan;
                % pd_unmerged{f}.data.Y(idxall)=nan;
                % pd_unmerged{f}.data.Z(idxall)=nan;
               %%%%%%%% 

                if ~exist("positionData","var")
                    positionData = pd_unmerged{f};
                else
                    positionData.data = [positionData.data; pd_unmerged{f}.data(:,{'X','Y','Z'})];
                    positionData.time = positionData.time + pd_unmerged{f}.time;
                end
            end


              %%%%%% normalize position %%%%%%%%
            positionData = positionData.normalizePositionData; % normalize to 100 by 100

            %%%%%% get angle for position %%%
            obj.positionAngle = positionData.getPositionDataAngle;

            %%% get turn times given angle %%%%%%%%%%%%%%%%
            %Turns.TurnLocation = location (in rad) when turn occurs
            %Turns.Index = index into positionData length vector
            %Turns.TimestampsZT = zt timestamp 
            %Turns.Direction = 0 if turns counterclockwise and 1 if turns clockwise
            obj.Events.Turns = position.directionsCircularTrackRadians(positionData);
          %  obj.Events.sensorData = obj.getSensorData;
          

        end

        function sensorData = getSensorData(obj)
            sensorpath = [obj.basepath '/SensorData']; %put here the txt files
            filenames = dir([sensorpath '/' '*.txt']);
            try
            for f = 1:length(filenames)
                
                fID = fopen([sensorpath '/' filenames(f).name],'r');
             %   T = fscanf(fID,['%11s'])


                fclose(fID)

                end
            catch
                disp('Sensor Data not available')
            end

        end

        %%%%%%%% to eliminate outside of physical track coordinates %%%%%%%
        function obj = filterPositionbyCoordinates(obj)

            positionData = positionData.eliminateOffCircTrack(circleParams); %CHECK with new sessions/animals

        end

        %%%% To use sleep scores to filter position %%%%%%%%%%%%%%%%%%%%%%%
        function obj = filterPositionbySleepScore(obj)

            % get states for entire LFP trace, as having sleep included
            % will make for better scoring

          sleepScore = callSleepScoreMaster(obj); %divide into rem,nrem, wake, and quiet wake
          disp('test')

          % get timestamps
          %sleep score is roughly 1 Hz but padded so to start I'll
          %guestimate 1Hz

          obj.LFP.TimeIntervalCombined.getTimePointsInAbsoluteTimes;
          sleepTimestamps = obj.LFP.TimeIntervalCombined.getTimeSeriesDownsampled(1250);

          lfpTime1Hz = obj.LFP.TimeIntervalCombined.getDownsampled(1250);
          sleepTimestamps = lfpTime1Hz.getTimePointsInAbsoluteTimes;% 3 points longer than the other

          % gonna chop the end for now
          sleepTimestamps = sleepTimestamps(1:end-3);
          posTimestamps = obj.positionData.time.getTimeStamps;
          
          positionTimesIdx = interp1(sleepTimestamps);
          for i = 1:length(posTimestamps)
              posIdx(i) = find(positionTimesIdx == posTimestamps(i));
          end
        end

        %%% Get sleep score %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function sleepScore = callSleepScoreMaster(obj)
        % Dividing into rem, nrem, quiet wake, active wake

        sleepScore = SleepScoreMaster(obj.basepath,'rejectchannels',obj.badChannels,...
            'EMGChannels',obj.EMGChannels,'Notch60Hz',1,'ThetaChannels',obj.thetaChannel,'overwrite',true);
        end

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