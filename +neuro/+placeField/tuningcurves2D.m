function [] = tuningcurves2D(varargin)

p = inputParser;
addRequired(p,'behaviortimestamps',@isstruct)
addRequired(p,'phyoutputpath',@ischar)
addParameter(p,'behaviorsamprate',60,@isnumeric)
addParameter(p,'optitrack',[],@isstruct) % NEED THIS OR TRIALS
addParameter(p,'trials',[],@istable) 
addParameter(p,'labelstouse',[1 2],@isnumeric)
addParameter(p,'fBase',1.0,@isnumeric) % The assumed basal firing rate for each unit
addParameter(p,'minT_occ',1.5,@isnumeric) % Minimum sampling time (in seconds) at each location
addParameter(p,'tau',0.25,@isnumeric) % The time-bin size (in seconds)
addParameter(p,'showplots',true,@islogical)

parse(p,varargin{:})
behaviortimestamps = p.Results.behaviortimestamps;
phyoutputpath = p.Results.phyoutputpath;
behaviorsamprate = p.Results.behaviorsamprate;
optitrack = p.Results.optitrack;
trials = p.Results.trials;
labelstouse = p.Results.labelstouse;
fBase = p.Results.fBase;
minT_occ = p.Results.minT_occ;
tau = p.Results.tau;
showplots = p.Results.showplots;

%% preprocessing

if isempty(optitrack) && isempty(trials)
    error('Need either optitrack or trials')
end

if ~isempty(optitrack) && ~isempty(trials)
    error('Can only have one of optitrack or trials')
end

t = behaviortimestamps.fromlfpfilestart/1000; % in sec

% right here, put in just the parts from the trials to have trial only 

if isempty(trials)
x = optitrack.position.interpolatedx;
y = optitrack.position.interpolatedy;
else
    % concatenate the trialdata
    x_alltrials = [];
    y_alltrials = [];
    for t = 1:size(trials,1)
        x_alltrials = horzcat(x_alltrials,trials.xtrialdata{t});
        y_alltrials = horzcat(y_alltrials,trials.ytrialdata{t});
    end
end

spikes = getspikes(phyoutputpath,'labelstouse',labelstouse); %,'lfptimestamps',lfptimestamps);

nUnits = length(spikes.cellIDs);
labels = unique(spikes.cellLabels);

%% Occupancy and tuning curves

[~,~,~,~,empties,FRmaps] = buildNBbasedTC02_RW(x,y,t,spikes,behaviorsamprate,...
    fBase,minT_occ,tau,showplots);

%% Plot 

if showplots == 1

    for iU = 1:nUnits

        if plot_on
            clf;
        end
        % figure
        spks = spikes.ms_spiketimes{iU}/1000; % into seconds
        xint = interp1(t,xb,spks,'nearest');   % Align the position of the animal when the cell spikes to an x-grid position
        yint = interp1(t,yb,spks,'nearest');   % Align the position of the animal when the cell spikes to an y-grid position

        spikeocc = zeros(nxbins,nybins);
        for iX = 1:nxbins
            % Find all X positions in the bin
            tmpX = (xint == iX);
            for iY = 1:nybins
                % Find all Y positions in the bin
                tmpY = (yint == iY);
                spikeocc(iX,iY) = sum(tmpX.*tmpY);  % Since tmpX and tmpY are binary, multiplying them will yield 1 only if both are 1.  This is effectively an AND function.
            end
        end
        smoothedspike = filter2(zz,spikeocc,'same');
        smoothedspike(empties) = nan;

        % prior, probability pre-data of being in the group
        % posterior, probability of being in the group given the data
        %
        alpha0 = fBase*minT_occ;    % alpha hyperparameter for gamma prior
        beta0 = minT_occ/tau;       % beta hyperparameter for gamma prior
        alpha = alpha0 + smoothedspike;   % alpha hyperparameter for gamma posterior
        beta = beta0 + smoothedocc/tau;   % beta hyperparameter for gamma posterior
        likelihood = alpha; % probability of spikes
        marg = beta./(beta + 1); % marginilization % prob of occupancy
        posterior = likelihood.*(1-marg)./marg; % p(A|B) = p(B|A).*(P(A)./p(B)    % called meanN
        %prob of alpha given beta, or prob of spikes given occupancy
        varN = likelihood.*(1-marg)./marg.^2;

        alphas(iU,:,:) = alpha;
        betas(iU,:,:) = beta;
        FRmap = smoothedspike./smoothedocc;     % This is your standard spikes per unit time based map.  It's not statistically correct.  But it's close and easy.
        FRmaps(iU,:,:) = FRmap;

        FRmap = FRmap';
        subplot(1,2,1)
        h = pcolor(FRmap);
        set(h,'linestyle','none')
        axis equal
        %         caxis([0 40]);
        colorbar
        axis off
        ans = find(FRmap == max(FRmap(:)));
        [xtmp ytmp] = ind2sub(size(FRmap),ans);
        %         hold on;
        %         plot(ytmp,xtmp,'wo','markersize',15);
        %         hold off;
        title('Raw mean FR (Hz)','fontsize',15);

        % Make spiking plot
        subplot(1,2,2);
        xint = interp1(t,x,spks,'nearest');   % Align the position of the animal when the cell spikes to an x-grid position
        yint = interp1(t,y,spks,'nearest');   % Align the position of the animal when the cell spikes to an y-grid position
        plot(x,y,'color',[0.5 0.5 0.5]);
        hold on;
        plot(xint,yint,'r.');
        title(['Unit: ' num2str(iU)],'fontsize',15);
        colorbar
        axis equal
        axis off
    end

end











