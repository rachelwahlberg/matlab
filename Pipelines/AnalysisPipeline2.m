classdef FigurePipeline

    %%FigurePipeline %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Rachel Wahlberg Aug 2023 %
    % This pipeline is largely based off of Utku Kaya code %

    properties
        basepath = '/data/ExperimentsRSW/CircularMaze/20230614/merged_20230614';
        basename = 'merged_20230614';
        behaviorpath = fullfile(basepath, 'behaviorfiles/');
        phyoutputpath =  fullfile(basepath,'home/wahlberg/merged_M1_20211124crs.GUI');
        TimeIntervalCombinedpath = fullfile(basepath, [basename '_1250Hz.TimeIntervalCombined.csv']);
        phyfolder = '/home/wahlberg/merged_20230614crs.GUI/';
        optifolder = '/data/ExperimentsRSW/CircularMaze/20230614/Optitrack';
    end

    methods(Access=private)

        function obj = FigurePipeline()

        end
    end

    methods(Static) % static means you can call it without first creating an instance of the class.
        % Concrete implementation.  See Singleton superclass.
        function obj = instance()
            persistent uniqueInstance %only one instance can be created at a time
            if isempty(uniqueInstance) % if no instance yet created
                obj = neuro.spike.SpikeFactory(); % create one
                uniqueInstance = obj;
            else
                obj = uniqueInstance;
            end
            %% 1a) Define the basepath of the dataset to run.
            % The dataset should at minimum consist of the raw data and spike sorted data.
            cd(basepath)
        end
    end

    methods

        %% METHODS
        function [] = callMethods
            filepathtosave = [basepath '/Figures/Methods/methods_' date '.fig'];
            createMethodsFigure(filepathtosave)
        end

        function [] = callOnlinePlaceCells
            %% ONLINE PLACE CELLS
            % a) ex place cell at beginning of task, end of pair 1, end of pair 2
            % b) heat map place cell as a function of time
            % c) fidelity?
            % d) jumpiness?
            % e) etc...
        end

        function [] = callOnlineTheta
            %% ONLINE THETA POWER
            % a) th power as func of time
            % b) speed as func of time
            % c) example locked raster over time
            % d) performance as func of time
            % e) summary stats...
            % f) more summary stats ...
            % g) look ahead distance
            % h) other redish comparisons?
        end

        function [] = callSWRintro
            %% SWRS + REPLAYS INTRO
            % a) SWR power vs #replays
            % b) example replay
            % c) ripple power vs performance
            % d) etc...
        end

        function [] = callReplayOverTime
            %% REPLAYS OVER TIME
            % a) examples over time
            % b) ask pho for inspo (with his long vs short track stuff)
        end

        function [] = callColorAnalyses
            %% COLOR ANALYSES
            % a) task rule (drawing)
            % b) example place cell at boundary
            % c) violin plots of #PC at boundary in cross versus same color pairs
            % d) black performance vs white vs cross
            % e) time to performance for b vs w vs c
            % f) etc...
        end

        function [] = callDentateGyrusAnalyses
            %% DG PRELIMINARY ANALYSES (second paper?)
            % a) example cells
            % b) violin plots summary stats
            % c) PSD
            % d) raster
            % e) etc... (go through hpc analyses but with dg cells)

        end

    end

end