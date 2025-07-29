%figureScript

%% 0a) merge files
% Run only pre-sorting! 
usedate = 20230614;
basefolders = preprocess.getbasefolders(usedate); % adapt to include the variables below to not have to rerun each time
mergedfilename = '/data/ExperimentsRSW/CircularMaze/Harry/CA1/20230614_2/merged_20230614';
preprocess.OpenEphysPreprocess_RW(basefolders,mergedfilename);

%% BEFORE SPIKE SORTING.
% get artifact periods.

mf = methodsFigure();
[cleanlfp,removetimestamps] = preprocess.glitchdetector(mf.LFP.Values',mf.LFP.TimeIntervalCombined.getTimePointsInSamples',...
    mf.LFP.TimeIntervalCombined.getSampleRate,'winsec',.1,'eventfilename','merged_20230614.bad.evt');

% Then spike sort

% then convert to res/clu files

Phy2Neurosuite(mf.basepath,mf.phypath)


%% Ripples
swr = neuro.ripple.SWRDetectionMethodSWR(mf.basepath);


% to RUN figurepipeline etc. 

mf.exampleTrace([3600 3605],[0.5 500])



as = AcrossSessions();


% timeWindow1 = [obj.positionData.time.getStartTime obj.positionData.time.getEndTime]; % for aligning to show same window

timeWindow = [obj.positionData.time.getStartTime + minutes(45) obj.positionData.time.getStartTime + minutes(50)];          
mf.testAlignment(timeWindow);


pcs = [1 2; 3 4];
%nonpcs = 2;
figure
hold on
bar(pcs,'stacked')

%Figure saving
h = gcf;
h.Renderer = 'painters';
figname = 'theta2_anesthesiology2023.pdf';
print(figname,'-dpdf','-r300','-bestfit'); %did smoosh the raster a bit tho





