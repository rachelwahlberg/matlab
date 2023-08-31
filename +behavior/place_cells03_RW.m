function [FRmaps,xbins,ybins,empties,alphas,betas] = place_cells03_RW(varargin)

%Adam Johnson
%Adapted Rachel Wahlberg 2022
p = inputParser;
addRequired(p,'behaviortimestamps',@isstruct)
addRequired(p,'optitrack',@isstruct)
addRequired(p,'phyoutputpath',@ischar) % folder for p
addParameter(p,'labelstouse',1:9,@isnumeric) % indexing into q in clusterinfo
addParameter(p,'basepath',[],@ischar)
addParameter(p,'behaviorsamprate',60,@isnumeric)
addParameter(p,'placecellplots',true,@islogical)
addParameter(p,'fBase',1.0,@isnumeric) % The assumed basal firing rate for each unit
addParameter(p,'minT_occ',1.5,@isnumeric) % Minimum sampling time (in seconds) at each location
addParameter(p,'tau',0.25,@isnumeric) % The time-bin size (in seconds)
addParameter(p,'savefigs',false,@islogical)

parse(p,varargin{:})
behaviortimestamps = p.Results.behaviortimestamps;
optitrack = p.Results.optitrack;
phyoutputpath = p.Results.phyoutputpath;
labelstouse = p.Results.labelstouse;
basepath = p.Results.basepath;
behaviorsamprate = p.Results.behaviorsamprate;
placecellplots = p.Results.placecellplots;
fBase = p.Results.fBase;
minT_occ = p.Results.minT_occ;
tau = p.Results.tau;
savefigs = p.Results.savefigs;

%% Preprocessing
% Get the (x,y) positions, their times t, the sampling rate for the video
% data, the number of units, and the average firing rate for each cell.
% [t,x,y,sample_rate,nUnits,aveFR] = loadTXY(unitdata);

t = behaviortimestamps.fromlfpfilestart/1000; % in sec
x = optitrack.position.interpolatedx;
y = optitrack.position.interpolatedy;

spikes = getspikes(phyoutputpath,'labelstouse',labelstouse); %,'lfptimestamps',lfptimestamps);

nUnits = length(spikes.cellIDs);
labels = unique(spikes.cellLabels);

%% Occupancy and tuning curves
% see definitions above
[alphas,betas,xbins,ybins,empties,FRmaps] = buildNBbasedTC02_RW(x,y,t,spikes,behaviorsamprate,...
    fBase,minT_occ,tau,placecellplots);

%% STOP HERE!!!!  We'll get to the other stuff later
%         alpha0 = fBase*minT_occ;    % alpha hyperparameter for gamma prior
%         beta0 = minT_occ/tau;       % beta hyperparameter for gamma prior
%         alpha = alpha0 + smoothedspike;   % alpha hyperparameter for gamma posterior
%         beta = beta0 + smoothedocc/tau;   % beta hyperparameter for gamma posterior

%% Plot rate maps
placecellplots = 0;
if placecellplots
    figure
    for iU = 1:nUnits
        alpha = squeeze(alphas(iU,:,:));
        beta = squeeze(betas(iU,:,:));
        r = alpha;
        p = beta./(beta + 1);
        meanN = r.*(1-p)./p;
        varN = r.*(1-p)./p.^2;
        subplot(2,1,1);
        h = pcolor(meanN'/tau);
        set(h,'linestyle','none')
        axis equal
        colorbar
        axis off
        title(['NB-based mean FR (Hz) - Unit:' num2str(iU)],'fontsize',15);
        colormap('jet');
        subplot(2,1,2);
        FRmap = squeeze(FRmaps(iU,:,:));
        h = pcolor(FRmap');
        set(h,'linestyle','none')
        axis equal
        colorbar
        axis off
        title('Poisson-based mean FR (Hz)','fontsize',15);
        colormap('jet');
        pause;
    end
end

% 
% %% Make non-spatial tuning curves
% %[output] = get_behavior01(samples);
% plot_on = 0;
% % tau = 0.5;
% % fBase = 1.5;
% % minT_occ = 1;
% [scounts,tcounts,levels,alphasNS,betasNS] = nonspatial_tc03(unitdata,samples,tau,fBase,minT_occ,plot_on);
% 
% 
% %% Run decoding analysis
% [output] = get_behavior01(samples);                 % Get behavior
% %%
% iS = 25;        % 11, 12, 15, 16, 18
% % Create the population vectors
% t_start = output(iS,1) - 10;                     % Set the start time for decoding
% t_stop = output(iS,1) + 1;              % Set the end time for decoding
% dt = tau;                           % This is the interval size again
% [popvec] = buildpopvec(unitdata,t_start,t_stop,dt,find(cells_use));
% % [popvec] = make_theta_popvec01(unitdata,theta_c);
% 
% [pSpace,xb,yb,t_use] = spatial_NBdecoding01(popvec,alphas,betas,x,y,t,t_start,t_stop,dt,xbins,ybins,empties,nUnits_use);
% [pNS] = nonspatial_NBdecoding01(popvec,alphasNS,betasNS,cells_use);
% 
% for iT = 1:size(popvec,1),
%     
%     
%     subplot(3,3,1:6);
%     h = pcolor(squeeze(pSpace(iT,:,:))');
%     hold on;
%     plot(xb(iT),yb(iT),'wo','markersize',10,'linewidth',2);
%     text(2,length(ybins) - 2,'1','fontsize',15,'color','k')
%     text(2,2,'3','fontsize',15,'color','k')
%     text(length(xbins) - 2,length(ybins) - 2,'2','fontsize',15,'color','k')
%     text(length(xbins) - 2,2,'4','fontsize',15,'color','k')
%     hold off;
%     set(h,'linestyle','none')
%     axis equal
%     colorbar
%     axis off
%     colormap('jet'); colorbar;
%     title(['Trial: ' num2str(iS) '/' num2str(output(iS,4)) ', time ' num2str(output(iS,1) - t_use(iT)) 's to choice']);
%     
%     subplot(3,3,7);
%     h = bar(pNS{1}(iT,:));
%     set(h,'facecolor',[0.5 0.5 0.5]);
%     ylabel('p(context)');
%     set(gca,'ylim',[0 1]);
%     title(['Context: ' num2str(output(iS,8))]);
%     subplot(3,3,8);
%     h = bar(pNS{2}(iT,:));
%     set(h,'facecolor',[0.5 0.5 0.5]);
%     ylabel('p(location)');
%     set(gca,'ylim',[0 1]);
%     title(['Location: ' num2str(output(iS,9))]);
%     subplot(3,3,9);
%     h = bar(pNS{3}(iT,:));
%     set(h,'facecolor',[0.5 0.5 0.5]);
%     ylabel('p(object)');
%     set(gca,'ylim',[0 1]);
%     title(['Object: ' num2str(output(iS,11)) '/Reward: ' num2str(output(iS,5))]);
%     pause%(0.05);
% end
% 
