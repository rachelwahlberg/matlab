function [spiketimes,spikeclusters,cells,cellCh,lab] = getspikeraster(varargin)

%% Get inputs
p = inputParser;
addRequired(p,'phyoutputpath',@ischar) %GUI folder
addRequired(p,'rawlfpsamprate',@isnumeric)
addParameter(p,'removebadtrials',false,@logical) % [starttime endtime] ms relative to start of
addParameter(p,'trials',[],@istable) % [starttime endtime] in ms
addParameter(p,'trialsPlottype','alltrials',@ischar)
addParameter(p,'plotbyfield',[],@ischar)
addParameter(p,'labelstouse',[1 2],@isnumeric) % what q labels to use - which cell types
addParameter(p,'showplots',true,@islogical)
addParameter(p,'newplotfolder',[],@ischar) % if included, will save figures
addParameter(p,'newplotbasename',[],@ischar) % need this as well to save figs

parse(p,varargin{:})
phyoutputpath = p.Results.phyoutputpath;
rawlfpsamprate = p.Results.rawlfpsamprate;
removebadtrials = p.Results.removebadtrials;
trials = p.Results.trials;
trialsPlottype = p.Results.trialsPlottype;
plotbyfield = p.Results.plotbyfield;
labelstouse = p.Results.labelstouse;
showplots = p.Results.showplots;
newplotfolder = p.Results.newplotfolder;
newplotbasename = p.Results.newplotbasename;

%% get data
% spikeclustersfilename = fullfile(basepath,[basename '_phy/'],'spike_clusters.npy');
% spiketimesfilename = fullfile(basepath,[basename '_phy/'],'spike_times.npy');
% clusterinfofilename = fullfile(basepath,[basename '_phy/'],'cluster_info.tsv');

raw_spiketimes = double(readNPY(fullfile(phyoutputpath,'spike_times.npy')));% NPY output from phy, spike_times
raw_spikeclusters = readNPY(fullfile(phyoutputpath,'spike_clusters.npy')); % NPY output from phy, spike_clusters
clusterinfo = tdfread(fullfile(phyoutputpath,'cluster_info.tsv')); % tsv output from phy, cluster_info

raw_mstimes = raw_spiketimes/rawlfpsamprate*1000;% get into ms

% cstart = 3729251;%*30000;
% cstop = 3730164;%*30000;
% inrippleSpikes = raw_mstimes>cstart & raw_mstimes<cstop;
% idx = find(inrippleSpikes == 1);

%% Divide by trial/type of trial
% Also make an option to only plot the clean trials
if isnumeric(clusterinfo.q(1))
    label = clusterinfo.q(:);
else
    label = zeros(length(clusterinfo.q),1);
    for clu = 1:length(clusterinfo.q)
        tmp = str2double(clusterinfo.q(clu,1));
        label(clu,1) = tmp; % choosing the label (for pyramidal, interneuron, mua, etc)
       % ch(clu,1) = clusterinfo.ch(clu);
    end
end

%select cells
c = 1;
for l = 1:length(labelstouse)
    idx = find(label == labelstouse(l));
    cells(c:c+length(idx)-1) = clusterinfo.id(idx);
    lab(c:c+length(idx)-1) = labelstouse(l);
    cellCh(c:c+length(idx)-1) = clusterinfo.ch(idx);
    c = c+length(idx);
end

if isempty(trials)
    spikeclusters{1} = raw_spikeclusters; % cell to make compatible with trial format
    spiketimes{1} = raw_mstimes;
    disp('summing across all trials')
else

    fname = plotbyfield; %field name to divide by
    if ~isempty(fname)
alltypes = trials.(fname);
plotby = unique(alltypes);
    end
end

switch trialsPlottype
    %%%%%%%%%%%%%%%% Plot all trials %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 'alltrials' % avg cells!
        startRaster{1} = trials.startTrial_file; % adding cell 1
        endRaster{1} = trials.endTrial_file;
        trialType = [];

        keep = 1;
        for tri = 1:length(startRaster{1})

            if removebadtrials == 1
                if trials.cleanTrial(tri) == 0
                    continue
                end
            end

            timeindex = raw_mstimes > startRaster{1}(tri) & ...
                raw_mstimes < endRaster{1}(tri);
            spiketimes{1}{keep} = raw_mstimes(timeindex);
            spikeclusters{1}{keep} = raw_spikeclusters(timeindex);
            keep = keep + 1;
        end

        ntypes = length(spiketimes);

    case 'allcells_bytrial' % plot a figure for each trial, with all cells
        %startRaster{1} = trials.startTrial_file; % adding cell 1
        %endRaster{1} = trials.endTrial_file;
        %trialType = [];

        [spiketimes,spikeclusters,ntypes] = get_trialspiketimes; %
        spiketimesArray = getspiketimesarray;

       % maxtrilength = max(cell2mat(cellfun(@length,spiketimes{1},'UniformOutput',false)));

        for type = 1:ntypes
            for tri = 1:length(spiketimes{type})
                for c = 1:length(cells)
                    cIdx = find(spikeclusters{type}{tri} == cells(c));
                    cSp{type}{tri}{c} = spiketimes{type}{tri}(cIdx);
                end
            end
        end

        % sort by FR

        len = cellfun('length',cSp{1}{5});
        [val,idx] = sort(len,'descend');

        figure
        hold on
        for type = 1:ntypes
            for tri = 1:length(spiketimes{type})
                title(['Spike Raster for Trial ' num2str(tri)])
                ylabel('Cell')
                xlabel('time (ms) from trial start')
                
                for c = 1:length(cSp{1}{1})
                    cIdx = idx(c);
                   %relativetime = cSp{type}{tri}{cIdx}/1000 - trials.startTrial_file(tri)/1250;
                     scatter(cSp{type}{tri}{cIdx}/1000,ones(length(cSp{type}{tri}{cIdx}),1)*c,"|","k",'LineWidth',1.2);
                   % scatter(relativetime,ones(length(cSp{type}{tri}{cIdx}),1)*c,"|","k",'LineWidth',1.2)
                end

                % just for trial 5!
                xlim([1 3100])
                pause
                clf
                hold on
            end
        end

% 
%         figure
%         hold on
%         for type = 1:ntypes
%             for tri = 1:length(spiketimes{type})
%                 title(['Spike Raster for Trial ' num2str(tri)])
%                 ylabel('Cell')
%                 xlabel('time (ms) from trial start')
%                 
%                 for c = 1:length(cSp{1}{1})
%                     cIdx = idx(c);
%                    % relativetime = cSp{type}{tri}{cIdx}/1000 - trials.startTrial_file(tri)/1250;
%                    
%                     scatter(relativetime,ones(length(cSp{type}{tri}{cIdx}),1)*c,"|","k",'LineWidth',1.2)
%                 end
% 
%                 % just for trial 5!
%                 xlim([1 3100])
%                 pause
%                 clf
%                 hold on
%             end
%         end

        %%%%%%%%%%%%%%%%%%%% for including in preliminary data fig
        %gcf
        %subplot(3,2,5:6)
        %for tri = 1:length(spiketimes{type})
        %tri = ;
%                 title(['trial ' num2str(tri)])
%                 ylabel('Cell')
%                 xlabel('t (ms) from trial start')
% 
%                 for c = 1:length(cSp{1}{1})
%                     relativetime = cSp{type}{tri}{c} - trials.startTrial_file(tri);
%                     scatter(relativetime,ones(length(cSp{type}{tri}{c}),1)*c,"|","k")
%                 end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%                
%             
%                 celltrials{c}{type} = logical(zeros(length(spiketimes{type}),maxtrilength));
% 
%                 
%                     celltrials{c}{type}(tri,1:size(spiketimesArray{type}{tri},2)) = logical(spiketimesArray{type}{tri}(c,:));
% 
%                     cell_trialsums{c}{type}(tri,1) = sum(celltrials{c}{type}(tri,:));
%                     cell_fr{c}{type}(tri,:) = cell_trialsums{c}{type}(tri,1)/length(spiketimesArray{type}{tri});
% 
%                 end
%                 cell_meanfr(c,type) = mean(cell_fr{c}{type});
% 
% 
%         if showplots == 1
% 
%             figure
% 
%             for c = 1:length(cells)
% 
%                 hold all
%                 ti = strcat('cell ', num2str(cells(c)));
%                 sgtitle(ti);
% 
%                 for type = 1:ntypes
%                     subplot(1,ntypes,type)
%                     plotSpikeRaster(celltrials{c}{type},'PlotType','vertline')
%                     %title([alltypes{type}])
%                 end
% 
%                 if ~isempty(newplotfolder)
%                     newplotname = fullfile(newplotfolder,[newplotbasename '.' num2str(cell(c)) '.fig']);
%                     save(gcf,newplotname);
%                     newplotname = fullfile(newplotfolder,[newplotbasename '.' num2str(cell(c)) '.pdf']);
%                     save(gcf,newplotname);
%                 end
%                 pause
%                 clf
%                 hold off
%             end
% 
%             figure
%             hold all
%             title('Mean FR per cell')
%             h = bar(cell_meanfr);
%             xlabel('cell')
%             ylabel('mean FR')
%             set(h,{'DisplayName'},{'Left Trials','Right Trials'}')
%             legend()
% 
%         end


            %%%%%%%%%%%%%%%% Sum Across Cells, split by trial condition %%%%%
        case 'condition_acrosscells'

            fname = plotbyfield; %field name to divide by
            alltypes = trials.(fname);
            plotby = unique(alltypes);

            %% Get relevant spiketimes
            % don't need to save the startTri/endTri if get the clusters at the same
            % time. Although maybe yes for the xaxis?

            [spiketimes,spikeclusters,ntypes] = get_trialspiketimes;
            spiketimesArray = getspiketimesarray;

            %% Plot raster
            %horzline (default) does not work! Needs to have the 'vertline' argument.
            %Otherwise plot is all blank, not skipping over nans.
            % which time format is this in?
            % v = 3

            for t = 1:ntypes
                ntrials = length(spiketimes{t});

                % for creating subplots
                if mod(ntrials,3) == 0; h = ntrials/3;
                elseif mod(ntrials,3) == 1; h = (ntrials+2)/3;
                elseif mod(ntrials,3) == 2; h = (ntrials+1)/3;
                end
                figure
                hold all
                ti = strcat(plotby(t),'trials');
                sgtitle(ti{1});
                for tri = 1:ntrials
                    subplot(3,h,tri)
                    plotSpikeRaster(spiketimesArray{t}{tri},'PlotType','vertline')
                end

                if ~isempty(newplotfolder)
                    newplotname = fullfile(newplotfolder,[newplotbasename '.' alltypes{t} 'trials.fig']);
                    save(gcf,newplotname);
                    newplotname = fullfile(newplotfolder,[newplotbasename '.' alltypes{t} 'trials.pdf']);
                    save(gcf,newplotname);
                end

                hold off
            end

            %%%%%%%%%%%%%%%% Individual cells, split by trial condition %%%%%
        case 'condition_withincells'

            [spiketimes,spikeclusters,ntypes] = get_trialspiketimes; %
            spiketimesArray = getspiketimesarray;

            maxtrilength_L = max(cell2mat(cellfun(@length,spiketimes{1},'UniformOutput',false)));
            maxtrilength_R = max(cell2mat(cellfun(@length,spiketimes{2},'UniformOutput',false)));
            maxtrilength = max(maxtrilength_L,maxtrilength_R); % max trial length altogether

            for c = 1:length(cells)
                for type = 1:ntypes
                    celltrials{c}{type} = logical(zeros(length(spiketimes{type}),maxtrilength));

                    for tri = 1:length(spiketimes{type})
                        celltrials{c}{type}(tri,1:size(spiketimesArray{type}{tri},2)) = logical(spiketimesArray{type}{tri}(c,:));

                        cell_trialsums{c}{type}(tri,1) = sum(celltrials{c}{type}(tri,:));
                        cell_fr{c}{type}(tri,:) = cell_trialsums{c}{type}(tri,1)/length(spiketimesArray{type}{tri});

                    end
                    cell_meanfr(c,type) = mean(cell_fr{c}{type});
                end
            end

            if showplots == 1
                figure
                for c = 1:length(cells)

                    hold all
                    ti = strcat('cell ', num2str(cells(c)));
                    sgtitle(ti);

                    for type = 1:ntypes
                        subplot(1,ntypes,type)
                        plotSpikeRaster(celltrials{c}{type},'PlotType','vertline')
                        title([alltypes{type}])
                    end

                    if ~isempty(newplotfolder)
                        newplotname = fullfile(newplotfolder,[newplotbasename '.' num2str(cell(c)) '.fig']);
                        save(gcf,newplotname);
                        newplotname = fullfile(newplotfolder,[newplotbasename '.' num2str(cell(c)) '.pdf']);
                        save(gcf,newplotname);
                    end
                    pause
                    clf
                    hold off
                end

                figure
                hold all
                title('Mean FR per cell')
                h = bar(cell_meanfr);
                xlabel('cell')
                ylabel('mean FR')
                set(h,{'DisplayName'},{'Left Trials','Right Trials'}')
                legend()



                %                 for c = 1:length(cells)
                %
                %                     hold all
                %                     ti = strcat('cell ', num2str(cells(c)));
                %                     sgtitle(ti);
                %
                %                     for type = 1:ntypes
                %                         subplot(1,ntypes,type)
                %                         histcounts(cell_meanfr)
                %                         title([alltypes{type}])
                %                     end
                %
                %                     if ~isempty(newplotfolder)
                %                         newplotname = fullfile(newplotfolder,[newplotbasename '.' num2str(cell(c)) '.fig']);
                %                         save(gcf,newplotname);
                %                         newplotname = fullfile(newplotfolder,[newplotbasename '.' num2str(cell(c)) '.pdf']);
                %                         save(gcf,newplotname);
                %                     end
                %                     pause
                %                     clf
                %                     hold off
                %                 end
            end
end
% cells([25,29,30]) = []; % temporarily - quite active

%%%%%%%%%%%%%%%%%%%%%%%%% get_trialspiketimes function %%%%%%%%%%%%%%%%%%%%

    function [spiketimes,spikeclusters,ntypes] = get_trialspiketimes()

        for t = 1:length(plotby)
            
            if ischar(plotby)
                type = plotby{t};
                typeindex = find(strcmp(alltypes(:),type) == 1);
                ntype = length(typeindex);
            elseif isnumeric(plotby)% didn't test this
                type = plotby(t);
                typeindex = find(alltypes(:) == type);
                ntype = length(typeindex);
            end

            keep = 1;
            for tri = 1:length(typeindex)

                if removebadtrials == 1
                    if trials.cleanTrial(tri) == 0
                        continue
                    end
                end

                startRaster{t}{keep} = trials.startTrial_file(typeindex(tri));
                endRaster{t}{keep} = trials.endTrial_file(typeindex(tri));

                timeindex = raw_mstimes > startRaster{t}{tri} & raw_mstimes < endRaster{t}{tri};
                spiketimes{t}{keep} = raw_mstimes(timeindex);
                spikeclusters{t}{keep} = raw_spikeclusters(timeindex);
                keep = keep + 1;
            end
        end % for t = 1:length(alltypes)

        ntypes = length(spiketimes);

    end

%%%%%%%%%%%%%%%%%%%%%%%%% cell_spiketimes_array function %%%%%%%%%%%%%%%%%%

    function spiketimesArray = getspiketimesarray() % logical 

        for type = 1:ntypes
            for tri = 1:length(spiketimes{type}) % if not including trials it will pass in the whole thing
                for c = 1:length(cells)
                    cellIDs{type}{tri}(c,:) = spikeclusters{type}{tri} == cells(c);
                    %     pyr_spiketimes{i,1} = (spike_times(spike_clusters==pyr(i)))';

                    cellIndex{type}{tri} = find(cellIDs{type}{tri}(c,:)==1);

                    trialSpiketimes{type}{tri}(c,:)=uint64(zeros(length(cellIDs{type}{tri}(c,:)),1));
                    trialSpiketimes{type}{tri}(c,cellIndex{type}{tri}) = ...
                        spiketimes{type}{tri}(cellIndex{type}{tri});
                    spiketimesArray{type}{tri}(c,:) = trialSpiketimes{type}{tri}(c,:) ~=0;
                    clear cellIndex
                end
            end
        end

    end

end

% check sampling rate
% check the overactive cells
% check the spike raster code overall






















