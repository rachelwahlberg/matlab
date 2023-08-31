function spikes = getspikes(varargin)
% pulls out spiketimes and spiking activity from phy outputs
% 1-3, pyramidal great,good,okay-ish
% 4-6, interneuron great,good,okay-ish
% 7-9, mua great,good,okay-ish

%spike_times.npy - [nSpikes, ] uint64 vector giving the spike time of each 
% spike in samples. To convert to seconds, divide by sample_rate from params.py.

% spike_clusters.npy - [nSpikes, ] int32 vector giving the cluster identity
% of each spike. This file is optional and if not provided will be automatically 
% created the first time you run the template gui, taking the same values as 
% spike_templates.npy until you do any merging or splitting

% cluster_info includes fields regarding cluster information - includes "q"
% which is self-created label column

p = inputParser;
addRequired(p,'phyoutputpath',@ischar)
addParameter(p,'rawlfpsamprate',30000,@isnumeric)
addParameter(p,'labelstouse',1:9,@isnumeric)
addParameter(p,'lfptimestamps',[],@isstruct)
parse(p,varargin{:})
phyoutputpath = p.Results.phyoutputpath;
rawlfpsamprate = p.Results.rawlfpsamprate;
labelstouse = p.Results.labelstouse;
lfptimestamps = p.Results.lfptimestamps;

%% get data
spikeclustersfile = fullfile(phyoutputpath,'spike_clusters.npy');
spiketimesfile = fullfile(phyoutputpath,'spike_times.npy');
clusterinfofile = fullfile(phyoutputpath,'cluster_info.tsv');

spiketimes = readNPY(spiketimesfile);
spikeclusters = readNPY(spikeclustersfile);
clusterinfo = tdfread(clusterinfofile);

% get into ms
%raw_mstimes = raw_spiketimes/rawlfpsamprate*1000;

%% get cell IDs for provided labels

if isnumeric(clusterinfo.q(1))
    label = clusterinfo.q(:);
else
    label = zeros(length(clusterinfo.q),1);
    for clu = 1:length(clusterinfo.q)
        tmp = str2double(clusterinfo.q(clu,1));
        label(clu,1) = tmp; % choosing the label (for pyramidal, interneuron, mua, etc)
    end
end

%select cells
c = 1;
for l = 1:length(labelstouse)
    ids = find(label == labelstouse(l));
    cells(c:c+length(ids)-1) = ids;
    lab(c:c+length(ids)-1) = labelstouse(l);
    c = c+length(ids);
end


% for clu = 1:length(clusterinfo.q)   
%     tmp = str2double(clusterinfo.q(clu,1));
%     label(clu,1) = tmp; % choosing the label (for pyramidal, interneuron, mua, etc)
% end
% 
% %select cells
% c = 1;
% for l = 1:length(labelstouse)
%     ids_ind = find(label == labelstouse(l));
%     ids = clusterinfo.id(ids_ind);
%     cell_ids(c:c+length(ids)-1) = ids;
%     lab(c:c+length(ids)-1) = labelstouse(l);
%     c = c+length(ids);
% end

%% Get spiking info for each cell
nsamps = length(spiketimes);
ncells = length(cells);
clustTimes = uint64(zeros(nsamps,ncells));

for c = 1:ncells
    clusters = spikeclusters == cells(c);
    clustTimes(:,c) = spiketimes.*uint64(clusters); % need?
    cellSpiketimes_samp{c} = spiketimes(clusters);
    cellSpiketimes_ms{c} = double(spiketimes(clusters))/rawlfpsamprate*1000; %di
end

%% Get average firing rate
% might be inaccurate due to merged files so commenting out for now
% nCells = length(cells);                    % The number of cells recorded
% avgFR = zeros(nCells,1);                            % Average firing rate - this is useful for distinguishing pyramidal cells from interneurons
% tstart = 0;
% tstop = length(spiketimes)/samplingrate; % in 
% for iU = 1:nCells,
%     avgFR(iU) = length(unitdata.units(iU).ts)/(tstop - tstart);     % This is the number of spikes divided by the total time
% end

%% Put info into structure

spikes.cellIDs = cells;
spikes.cellLabels = lab;
spikes.ms_spiketimes = cellSpiketimes_ms;
spikes.samp_spiketimes = cellSpiketimes_samp;
spikes.samplerate = rawlfpsamprate;







