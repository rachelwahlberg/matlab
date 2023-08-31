%folder='/home/wahlberg/Exp_Data/prelimdata_rat2_090221_1/prelimdata_rat2_090221_1/prelimdata_rat2_090221_1crs-1.GUI/';
 folder='/home/wahlberg/Exp_Data/prelimdata_rat2_090121_1/prelimdata_rat2_090121_1/prelimdata_rat2_090121_1crs.GUI/';

sf=neuro.spike.SpikeFactory.instance;
sa=sf.getPhyOutputFolder(folder);
sa=sa.sort('ch');
%%%% For specific timeframe %%%%%%
windowInSeconds=[100 101];
winInSamples=windowInSeconds*sa.TimeIntervalCombined.getSampleRate;
idx=sa.SpikeTable.SpikeTimes>winInSamples(1)&sa.SpikeTable.SpikeTimes<winInSamples(2);
sa1=sa;
sa1.SpikeTable=sa.SpikeTable(idx,:);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% sa=sa.getTimeInterval(1000,2000);
% try close(1);catch, end;figure(1);
sa1.plotRaster;