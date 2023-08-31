%%Process data - spike raster
clear
%get data
spike_clusters = readNPY('/home/wahlberg/Exp_Data/prelimdata_rat2_090221_1/prelimdata_rat2_090221_1/prelimdata_rat2_090221_1crs-1.GUI/spike_clusters.npy');
spike_times = readNPY('/home/wahlberg/Exp_Data/prelimdata_rat2_090221_1/prelimdata_rat2_090221_1/prelimdata_rat2_090221_1crs-1.GUI/spike_times.npy');
cluster_info = tdfread('/home/wahlberg/Exp_Data/prelimdata_rat2_090221_1/prelimdata_rat2_090221_1/prelimdata_rat2_090221_1crs-1.GUI/cluster_info.tsv','tab');

%get labels, select pyramidals level 1/2
label = str2num(cluster_info.q(:,1)); % 1 = best pyramidal, 2 = less good pyr, 3 = less good pyr, 6 = mua, 8 = interneuron
pyr = find(label == 1 | label == 2); % gives cluster numbers for the okay pyramidal neurons

for i = 1:length(pyr)
    pyr_clusters(i,:) = spike_clusters == pyr(i);
%     pyr_spiketimes{i,1} = (spike_times(spike_clusters==pyr(i)))';

    cluster_times = find(pyr_clusters(i,:)==1);

    pyr_spiketimes2(i,:)=uint64(zeros(length(pyr_clusters(i,:)),1));
    pyr_spiketimes2(i,cluster_times) = spike_times(cluster_times);
    pyr_spiketimes_array(i,:) = pyr_spiketimes2(i,:) ~=0;
    clear cluster_times
end

plotSpikeRaster(pyr_spiketimes_array) %line 196 - the plot function is not skipping over the nans in the data, it's making the 
%whole figure blank. There's for sure data in the structures and for
%example plot(xPoints(1:2),yPoints(1:2)) will plot but
%plot(xPoints(1:5),yPoints(1:5)) will not plot because of the nan. Not sure
%the issue? 

