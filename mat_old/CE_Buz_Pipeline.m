%% Cell Explorer/Buzcode Pipeline
% taking the best from them both

%  CE 1.Define the basepath of the dataset to run. The dataset should at minimum consist of the raw data and spike sorted data.
basepath = '/home/wahlberg/Analysis/M1_Nov2021/20211123/merged_M1_20211123/';
basename = 'merged_M1_20211123';
cd(basepath)

%% CE 2. Generate session metadata struct using the template function and display the meta data in a gui
session = sessionTemplate(pwd,'showGUI',false);

% You can also inspect the session struct by calling the GUI directly:
% session = gui_session(session);

% And validate the required and optional fields
validateSessionStruct(session);

%% BZ Get the lfp data from the dat file
% lfp is a struct
if isfile([basepath basename '_lfp.mat'])
    disp('Loading existing LFP file')
    load([basepath basename '_lfp.mat']);
else
    disp('LFP file does not yet exist, creating one from bz_GetLFP')
    bz_LFPfromDat(basepath)
    lfp = bz_GetLFP('all','basepath',basepath,'downsample',24); % or use m =memmapfile(merged_M1_20211123_raw.lfp','Format','int16') and then do the channels/division by SR manually
end

%% Remove artifacts
% RW function (in Code/mat)
[cleanlfp,removetimestamps] = glitchdetector(lfp.data,'eventfilename',[basepath basename, '.bad.evt']);






%% BZ Get SWR
good_lfp.data = double(lfp.data(3800000:4580000,[2:121,123:end])); % wrong
good_lfp.channels = [2:121,123:128];
ripchan = bz_GetBestRippleChan(good_lfp); %doesn't work if you have extremely noisy segments/bad channels included, so just select a good region of lfp

% Params (defaults listed below from buzcode tutorial
ripthresh = [2 5]; %min std, max std
ripdur = [30 400];%100 was default, 400 KD's suggestion
ripfreq = [130 200]; %200 was default, 250 KD's suggestion

% Detect SPW-Rs or load structure with detection info if already detected
ripples = bz_FindRipples(basepath,ripchan,'thresholds',ripthresh,...
    'durations',ripdur,'show','on','saveMat',true);

SaveRippleEvents([basename 'merged_M1_20211123_raw.evt.rip'],ripples,ripchan) %ripchan is 127
% Calculate ripple stats
[maps,data,stats] = bz_RippleStats(double(lfp.data),lfp.timestamps,ripples);

% Sort by duration vs amplitude of ripple
[~,dursort]=sort(data.duration,1,'descend');
[~,ampsort]=sort(data.peakAmplitude,1,'descend');

%% BZ 3. Intro spectral processing: filter in ripple frequency and compare to output of bz_FindRipples

% Filter channel in ripple freq range
    ripfiltlfp = bz_Filter(lfp.data,'passband',ripfreq,'filter','fir1'); % lots of options in bz_Filter, read through it
        
%% BZ 4.   Plots 
%% BZ 4i.  Look: Plot some ripples
    figure()
    x=maps.ripples(dursort,:);
    for i=1:100
        subplot(10,10,i)
        plot(x(i,:))
        set(gca,'XTick',[]);
        set(gca,'YTick',[]);
        axis off
    end

bad_timestamps = [28967,669060,2055840,3431210;127668,719530,2070740,3453060]

%% Get sleep states
bad_timestamps = [28967,127668;669060,719530;2055840,2070740;3431210,3453060];
bad_times = bad_timestamps/1000; % put into seconds
SleepState = SleepScoreMaster(basepath,'ignoretime',bad_times,'rejectchannels',[0 122]);
%SleepState = SleepScoreMaster_RW(basepath,'ignoretime',bad_timestamps,'rejectchannels',[0 122],'SWChannels',ripchan);
%% CE 3.1.1 Run the cell metrics pipeline 'ProcessCellMetrics' using the session struct as input
cell_metrics = ProcessCellMetrics('session', session,'showGUI',false);

%OR %%
cell_metrics = load("merged_M1_20211123_raw.dat.cell_metrics.cellinfo.mat");

%% CE 3.1.2 Visualize the cell metrics in CellExplorer
cell_metrics = CellExplorer('metrics',cell_metrics); 

%% CE 3.2 Open several session from basepaths
basepaths = {'/your/data/path/basename_1/','/your/data/path/basename_2/'};
cell_metrics = loadCellMetricsBatch('basepaths',basepaths);
cell_metrics = CellExplorer('metrics',cell_metrics);

%% Get place fields
basepath = '/home/wahlberg/Exp_Data/M1_Nov2021/20211123/merged_M1_20211123_raw'; 
clustering_path = '/home/wahlberg/Exp_Data/M1_Nov2021/20211123/merged_M1_20211123/merged_M1_20211123/merged_M1_20211123crs.GUI';
Phy2Neurosuite(basepath,clustering_path)
spikes = bz_GetSpikes('basepath',basepath,'saveMat',true);

% using optitrack format, but using the timestamps from behavior_timestamps
% (EDIT FOR CLARITY)
firingRateMap = generate_FiringRateMap_1D(optitrack2,spikes)





% %% 4. load a subset of units fullfilling multiple of criterium
% 
% % Get cells that are assigned as 'Interneuron'
% cell_metrics_idxs1 = loadCellMetrics('cell_metrics',cell_metrics,'putativeCellType',{'Interneuron'});
% 
% % Get cells that are have a tag 'InverseSpike' or 'Good' and are assigned as 'Interneuron'
% cell_metrics_idxs2 = loadCellMetrics('cell_metrics',cell_metrics,'tags',{'InverseSpike','Good'},'putativeCellType',{'Interneuron'});
%  
% 
% ## Flowcharts
% The flowcharts below show the processes in details. The boxes are color coded according to external files (grey), data structures (yellow), CellExplorer functions (green), and database (purple).
% 
% ### Preparing experimental metadata
% First step is gathering metadata in the session struct. This struct can contain all metadata necessary for calculating the cell metrics. You can use the `sessionTemplate` to extract and define the parameters and visualize them with the graphical interface `gui_session` for further manual entry. The templates will scan the basepath for specific files, and import existing metadata to minimize the manual entry. You can customize the template script to fit and extract information relevant to your data. [The session struct is defined here]({{"/datastructure/data-structure-and-format/#session-metadata"|absolute_url}}). The session struct follows the database structure of the Buzsaki Lab and all metadata can be loaded directly from the database for database sessions. See the example code below on how perform the actions in Matlab.
% ![](https://buzsakilab.com/wp/wp-content/uploads/2020/05/Flowcharts_Metadata.png){: .mt-4}
% 
% ### Processing cell_metrics
% Following the definition of metadata, the cell metrics calculation process can be performed. A single script processes all default cell_metrics (which can be customized and expanded). The process is fully automatic, except for the detection of monosynaptic connections, in which a graphical interface is shown for manual curation (the manual step can be turned off). See the [full list of default cell_metrics here]({{"/datastructure/standard-cell-metrics/"|absolute_url}}). Below follows two flowcharts: a simple with the minimal inputs and an advanced flowchart. The advanced chart shows all relevant files that are compatible, auto-detected and loaded by the cell_metrics calculation process.
% ![](https://buzsakilab.com/wp/wp-content/uploads/2020/05/Flowcharts_ProcessingModule.png){: .mt-4}
% 
% ### Running CellExplorer
% CellExplorer can be used with single recording sessions as well as batches of sessions. Batch loading is performed with the script `loadCellMetricsBatch`. The advanced flowchart below further details the capabilities of loading various GUIs from CellExplorer (`gui_session`, `gui_MonoSyn` and `gui_DeelSuperficial`) as well as do spike raster plots, that requires access to the local `spikes` struct and potentially also manipulation and events files when plotting PSTHs.
% ![](https://buzsakilab.com/wp/wp-content/uploads/2020/05/Flowcharts_GraphicalInterface.png){: .mt-4}
% 
% ## Running pipeline from a data path
% The pipeline follows the data standards [described here]({{"/datastructure/data-structure-and-format/"|absolute_url}}). Saving your data in the specified data formats, integrates your data better with CellExplorer, allowing you to plot spike rasters and event histograms among other things.
% 
% To run the pipeline from a session struct, please see this example
% [sessionTemplate.m](https://github.com/petersenpeter/CellExplorer/blob/master/calc_CellMetrics/sessionTemplate.m) file for how to format this properly. Please edit the template file to fit it to your data.
% ```m
% session = sessionTemplate;
% ```
% You can also view the session struct in a GUI:
% ```m
% session = gui_session(session);
% ```
% 
% To run the processing script from the Matlab Command Window from the session struct type:
% ```m
% cell_metrics = ProcessCellMetrics('session', session);
% ```
% You can also run it directly from a basepath and generate the session struct directly:
% ```m
% cell_metrics = ProcessCellMetrics;
% ```
% When calling the processing script with the sessionTemplate, a GUI will be shown allowing you to edit  metadata, both input parameters and the session struct. 
% 
% Once complete, view the result in CellExplorer by typing:
% ```m
% cell_metrics = CellExplorer('metrics',cell_metrics);
% ```
% ### Running CellExplorer in batch mode from list of data paths
% To open multiple sessions together you can run CellExplorer in batch mode. Below is an example for running CellExplorer on three sessions from the database:
% 
% ```m
% bsasepaths = {'sessionName1','sessionName2','sessionName3'};
% cell_metrics = loadCellMetricsBatch('basepaths',bsasepaths);
% cell_metrics = CellExplorer('metrics',cell_metrics);
% ```
% As you perform classifications in CellExplorer, you may save back to the original cell metrics stored with the sessions defined above. You can perform the batch mode from a list of paths as well.
% 
% ## Running pipeline using the Buzsaki lab database (for Buzsaki lab members)
% CellExplorer processing module `ProcessCellMetrics` uses a single Matlab struct for handling metadata. The struct is automatically loaded from the -buzsaki lab database if you are running the pipeline with the database, and is located in the base path once a session has been processed. To run the pipeline on a session named 'PetersSession' using the database type:
% ```m
% cell_metrics = ProcessCellMetrics('sessionName','PetersSession');
% ```
% To view the result in CellExplorer type:
% ```m
% cell_metrics = CellExplorer('metrics',cell_metrics);
% ```
% 
% ### Running CellExplorer in batch mode from database
% To open multiple sessions together you can run CellExplorer in batch mode. Below is an example for running CellExplorer on three sessions from the database:
% 
% ```m
% sessionNames = {'sessionName1','sessionName2','sessionName3'};
% cell_metrics = loadCellMetricsBatch('sessions',sessionNames);
% cell_metrics = CellExplorer('metrics',cell_metrics);
% ```
% As you perform classifications in CellExplorer in batch mode, you can save your progress to the original sessions. You can work in batch mode from a list of paths as well.
% 
% 
