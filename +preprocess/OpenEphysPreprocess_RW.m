function [] = OpenEphysPreprocess_RW(basefolders,mergedfilename)

%Utku Kaya, original directory:
%/home/wahlberg/Toolboxes/ephys/Analysis/MATLAB/Ephys

% INPUTS
% basefolders should be structure.oebin files for all recordings
% mergefilename is the name to save the new file to

% OUTPUTS
% basename.dat; _1250Hz.lfp; _1250Hz.Probe.xlsx; 
% .TimeIntervalCombined.csv; _TimeIntervalCombined.csv; Probe.csv

oef=openEphys.OpenEphysRecordFactory;
for ifolder=1:numel(basefolders)
    basefolder=basefolders{ifolder};
    oer=oef.getOpenEphysRecord(basefolder);
    if ifolder==1
        oerc=oer;
    else
        oerc=oerc+oer;
    end

    
end


dfcl=oerc.mergeBlocksOfChannels(1:128,mergedfilename);
ctdh=neuro.basic.ChannelTimeDataHard(dfcl.getDataFile);
ctdds=ctdh.getDownSampled(1250);
% ctdds.getChannel(20).plot;

