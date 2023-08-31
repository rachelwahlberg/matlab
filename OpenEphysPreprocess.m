function [] = OpenEphysPreprocess_RW(basefolders)

%Utku Kaya, original directory:
%/home/wahlberg/Toolboxes/ephys/Analysis/MATLAB/Ephys

% basefolders should be structure.oebin files for all recordings

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
dfcl=oerc.mergeBlocksOfChannels(1:128,'/home/wahlberg/Exp_Data/M1_Nov2021/20211123/merged_M1_20211123_raw');
ctdh=neuro.basic.ChannelTimeDataHard(dfcl.getDataFile);
ctdds=ctdh.getDownSampled(1250);
% ctdds.getChannel(20).plot;

