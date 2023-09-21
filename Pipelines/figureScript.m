%figureScript

%% 0a) merge files
% Run only pre-sorting! 
usedate = 20230615;
basefolders = preprocess.getbasefolders(usedate); % adapt to include the variables below to not have to rerun each time
mergedfilename = '/data/ExperimentsRSW/CircularMaze/Harry/CA1/20230615/merged_20230615';
preprocess.OpenEphysPreprocess_RW(basefolders,mergedfilename);




% to RUN figurepipeline etc. 

mf = methodsFigure();

mf.exampleTrace([3600 3605],[0.5 500])



as = AcrossSessions();





