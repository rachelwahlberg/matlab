function basefolders = getbasefolders(usedate)

%%%% 11/23/21 %%%%
if usedate == 112321
basefolders={'/home/wahlberg/Exp_Data/M1_Nov2021/20211123/2021-11-23_15-34-56/RecordNode108/experiment1/recording1/structure.oebin',...
    '/home/wahlberg/Exp_Data/M1_Nov2021/20211123/2021-11-23_15-34-56/RecordNode108/experiment2/recording1/structure.oebin',...
    '/home/wahlberg/Exp_Data/M1_Nov2021/20211123/2021-11-23_15-51-17/RecordNode108/experiment1/recording1/structure.oebin',...
    '/home/wahlberg/Exp_Data/M1_Nov2021/20211123/2021-11-23_16-00-39/RecordNode108/experiment1/recording1/structure.oebin',...
    '/home/wahlberg/Exp_Data/M1_Nov2021/20211123/2021-11-23_16-20-34/RecordNode108/experiment1/recording1/structure.oebin',...
    '/home/wahlberg/Exp_Data/M1_Nov2021/20211123/2021-11-23_16-41-47/RecordNode108/experiment1/recording1/structure.oebin',...
    };

%%%% 11/24/21 %%%%
elseif usedate == 112421
basefolders={'/data/20211124/2021-11-24_08-52-16/RecordNode108/experiment1/recording1/structure.oebin',...
    '/data/20211124/2021-11-24_11-12-47/RecordNode108/experiment1/recording1/structure.oebin',...
    '/data/20211124/2021-11-24_11-12-47/RecordNode108/experiment1/recording2/structure.oebin',...
    '/data/20211124/2021-11-24_11-23-03/RecordNode108/experiment1/recording1/structure.oebin',...
    '/data/20211124/2021-11-24_11-34-01/RecordNode108/experiment1/recording1/structure.oebin',...
    '/data/20211124/2021-11-24_11-34-01/RecordNode108/experiment1/recording2/structure.oebin',...
    '/data/20211124/2021-11-24_12-55-19/RecordNode108/experiment1/recording1/structure.oebin',...
    };
% elseif usedate == 112421
% basefolders={'/home/wahlberg/Analysis/M1_Nov2021/20211124/2021-11-24_08-52-16/RecordNode108/experiment1/recording1/structure.oebin',...
%     '/home/wahlberg/Analysis/M1_Nov2021/20211124/2021-11-24_11-12-47/RecordNode108/experiment1/recording1/structure.oebin',...
%     '/home/wahlberg/Analysis/M1_Nov2021/20211124/2021-11-24_11-12-47/RecordNode108/experiment1/recording2/structure.oebin',...
%     '/home/wahlberg/Analysis/M1_Nov2021/20211124/2021-11-24_11-23-03/RecordNode108/experiment1/recording1/structure.oebin',...
%     '/home/wahlberg/Analysis/M1_Nov2021/20211124/2021-11-24_11-34-01/RecordNode108/experiment1/recording1/structure.oebin',...
%     '/home/wahlberg/Analysis/M1_Nov2021/20211124/2021-11-24_11-34-01/RecordNode108/experiment1/recording2/structure.oebin',...
%     '/home/wahlberg/Analysis/M1_Nov2021/20211124/2021-11-24_12-55-19/RecordNode108/experiment1/recording1/structure.oebin',...
%     };

elseif usedate == 112521
%%%% 11/25/21 %%%%
basefolders={'/home/wahlberg/Exp_Data/M1_Nov2021/2021-11-25/2021-11-25_08-56-19/RecordNode108/experiment1/recording1/structure.oebin',...
             '/home/wahlberg/Exp_Data/M1_Nov2021/2021-11-25/2021-11-25_10-58-43/RecordNode108/experiment1/recording1/structure.oebin',...
             };
elseif usedate == 20230613
    basefolders ={'/data/ExperimentsRSW/CircularMaze/Harry/CA1/20230613/2023-06-13_15-03-08/RecordNode101/experiment1/recording2/structure.oebin',...
        '/data/ExperimentsRSW/CircularMaze/Harry/CA1/20230613/2023-06-13_15-03-08/RecordNode101/experiment1/recording2/structure.oebin',...
        '/data/ExperimentsRSW/CircularMaze/Harry/CA1/20230613/2023-06-13_15-03-08/RecordNode101/experiment1/recording3'};
elseif usedate == 061423 %OLD
basefolders={'/data/20230614/2023-06-14/2023-06-14_11-15-49/RecordNode101/experiment1/recording1/structure.oebin',...
'/data/20230614/2023-06-14/2023-06-14_12-02-19/RecordNode101/experiment1/recording1/structure.oebin',...
'/data/20230614/2023-06-14/2023-06-14_12-23-27/RecordNode101/experiment1/recording1/structure.oebin',...
'/data/20230614/2023-06-14/2023-06-14_13-08-38/RecordNode101/experiment1/recording1/structure.oebin'};
elseif usedate == 20230615
    basefolders={'/data/ExperimentsRSW/CircularMaze/Harry/CA1/20230615/2023-06-15/2023-06-15_10-17-51/RecordNode101/experiment1/recording2/structure.oebin',...
        '/data/ExperimentsRSW/CircularMaze/Harry/CA1/20230615/2023-06-15/2023-06-15_10-17-51/RecordNode101/experiment1/recording2/structure.oebin',...
        '/data/ExperimentsRSW/CircularMaze/Harry/CA1/20230615/2023-06-15/2023-06-15_10-17-51/RecordNode101/experiment1/recording2/structure.oebin',...
        '/data/ExperimentsRSW/CircularMaze/Harry/CA1/20230615/2023-06-15/2023-06-15_10-17-51/RecordNode101/experiment1/recording2/structure.oebin',...
        '/data/ExperimentsRSW/CircularMaze/Harry/CA1/20230615/2023-06-15/2023-06-15_10-17-51/RecordNode101/experiment1/recording2/structure.oebin',...
        '/data/ExperimentsRSW/CircularMaze/Harry/CA1/20230615/2023-06-15/2023-06-15_10-17-51/RecordNode101/experiment1/recording2.structure.oebin'};
else
    disp('manually add oebin filepaths to this function for the date selected')
end

% %% 1a) Define the basepath of the dataset to run. 
% % The dataset should at minimum consist of the raw data and spike sorted data.
% basepath = '/home/wahlberg/Analysis/M1_Nov2021/20211123/merged_M1_20211123';
% basename = 'merged_M1_20211123';
% behaviorpath = fullfile(basepath, 'behaviorfiles/');
% manualbehaviorpath = '/home/wahlberg/Exp_Data/M1_Nov2021/20211123/Behavioral/221123_manualscoring.csv';
% phyoutputpath =  fullfile(basepath,'merged_M1_20211123_phy');
% % NEED THIS OR
% baseLFPfolders={'/home/wahlberg/Analysis/M1_Nov2021/20211123/2021-11-23_15-51-17/RecordNode108',...
%     '/home/wahlberg/Analysis/M1_Nov2021/20211123/2021-11-23_16-00-39/RecordNode108',...
%     '/home/wahlberg/Analysis/M1_Nov2021/20211123/2021-11-23_16-20-34/RecordNode108',...
%     '/home/wahlberg/Analysis/M1_Nov2021/20211123/2021-11-23_16-41-47/RecordNode108'};
% % THIS
% TimeIntervalCombinedpath = fullfile(basepath, [basename '_1250Hz.TimeIntervalCombined.csv']);
