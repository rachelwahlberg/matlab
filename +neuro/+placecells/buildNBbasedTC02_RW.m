function [alphas,betas,xbins,ybins,empties,FRmaps] = buildNBbasedTC02_RW(x,y,t,spikes,behsamprate,...
    fBase,minT_occ,tau,plot_on)

% fBase = 1.0;                % The assumed basal firing rate for each unit
% minT_occ = 1.5;             % Minimum sampling time (in seconds) at each location
% tau = 0.25;                 % The time-bin size (in seconds)
if nargin < 6
    plot_on = 0;
end
nUnits = length(spikes.cellIDs);                    % The number of cells recorded

%% Occupancy and tuning curves
binsize = 2; %2 %20; % I believe in cm
xbins = min(x):binsize:max(x);      % The x-coordinate bins for the grid
ybins = min(y):binsize:max(y);      % The y-coordinate bins for the grid
xbins = xbins(1:end-1) + 0.5*binsize;
ybins = ybins(1:end-1) + 0.5*binsize;
nxbins = length(xbins);
nybins = length(ybins);
xbinInd = 1:nxbins;
ybinInd = 1:nybins;

% Find the occupancy - how long does the animal spend in each location?
occ = zeros(nxbins,nybins); 
xb = interp1(xbins,xbinInd,x,'nearest');
yb = interp1(ybins,ybinInd,y,'nearest');
for iX = 1:nxbins
    % Find all X positions in the bin
    xtmp = (xb == iX);
    for iY = 1:nybins
        % Find all Y positions in the bin
        ytmp = (yb == iY);
        occ(iX,iY) = behsamprate*sum(xtmp.*ytmp);
    end
end
empties = find(occ < 1);

%% Smoothing kernel 
%   We blur the position of the animal because we don't exactly where it
%   is.  The smoothing kernal basically spreads the position of the animal
%   at each video tracker time as a 2D Gaussian (normal) distribution with
%   a center where we believe the animal's head was according to the video.
%   4.2 pixel/cm
rr = 1;
[xx,yy] = meshgrid(-2*rr:2*rr,-2*rr:2*rr);
zz = exp(-((xx).^2 + (yy).^2)/4); %the filter
zz = zz/(sum(zz(:)));
% figure; surfl(xx,yy,zz)
smoothedocc = filter2(zz,occ,'same');
smoothedocc(empties) = nan;

if plot_on
    figure;
    h = pcolor(smoothedocc');
    axis equal
    colorbar
    title('Time spent in each position')
    set(h,'linestyle','none')
    axis equal
end

%% Display tuning curves - press any key to get to the next plot
alphas = nan(nUnits,size(occ,1),size(occ,2));
betas = nan(nUnits,size(occ,1),size(occ,2));
FRmaps = nan(nUnits,size(occ,1),size(occ,2));

if plot_on
    figure;
end
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
        
        if plot_on
            % Make average firing rate
            FRmap = FRmap';
            subplot(2,2,1)
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
            subplot(2,2,2);
            xint = interp1(t,x,spks,'nearest');   % Align the position of the animal when the cell spikes to an x-grid position
            yint = interp1(t,y,spks,'nearest');   % Align the position of the animal when the cell spikes to an y-grid position
            plot(x,y,'color',[0.5 0.5 0.5]);
            hold on;
            plot(xint,yint,'r.');
            title(['Unit: ' num2str(iU)],'fontsize',15);
            colorbar
            axis equal 
            axis off

            subplot(2,2,3)
            hold on
             plot(t,xb);
            if ~isempty(xint); scatter(1:length(xint),xint,'r','LineWidth',2); end;
            ylabel('X pos')
            xlabel('time')
            title('X versus time')
            
            subplot(2,2,4)
            hold on
            plot(t,yb);
            if ~isempty(yint); scatter(1:length(yint),yint,'r','LineWidth',2); end;
            ylabel('Y pos')
            xlabel('time')
            title('Y versus time')
           

%             
%             % Make proper posterior plot for
%             subplot(2,2,3)
%             h = pcolor(posterior'/tau);
%             set(h,'linestyle','none')
%             axis equal
%             colorbar
%             axis off
%             title('NB-based mean FR (Hz)','fontsize',15);
%             
%             subplot(2,2,4)
%             h = pcolor(FRmap - (posterior'/tau));
%             %         h = pcolor(sqrt(varN)'/tau);
%             set(h,'linestyle','none')
%             axis equal
%             colorbar
%             axis off
%             title('Raw Mean FR - NB-based mean FR (Hz)','fontsize',15);
%             colormap('jet')
            
%             set(gcf,'position',[34 1 1297 797])
            pause;
    end
end


