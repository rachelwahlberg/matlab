function [] = plotPlaceFields(plottype,plottitle,varargin)
% meant to make creating different combinations of place field plots
% easier.

%%%%% EXAMPLE USE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% in external function ...
%figure
%grid_height = 1; grid_width = 2;
%h1 = tiledlayout(grid_height,grid_width);
%position = 1; h = 1; w = 1; ax1 = nexttile(position,[h,w]);
%position = 2; h = 1; w = 1; ax2 = nexttile(position,[h,w]);

%axes(ax1)
%ttl = 'Spk Locations';
%plot.plotPlaceFields('spkLocs',ttl,x,y,xspikes,yspikes)

%axes(ax2)
%ttl = 'Smooth Occ. Map';
%plot.plotPlaceFields('occMap',ttl,timespent)

%%% SEE METHODSFIGURE THOUGH. Cause to get occ map, for example, you
%%% need to run PFClassic first. See MethodsFigure for example.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

switch plottype
    case 'spkLocations' %plot spike locations on top of all occupancy locations.
        x = varargin{1};
        z = varargin{2};
        xspikes = varargin{3};
        zspikes = varargin{4};

        hold on
        title(plottitle)
        plot(x,z,'y')
        plot(xspikes,zspikes,'.r')
        ax = gca;
        ax.DataAspectRatio = [1 1 1];
        hold off

    case 'imagescMap' %applies to occupancy map, spikemap, etc
        data = varargin{1}; %TimeSpent, nSpikes, etc
        
        title(plottitle)
        hold on
        imagesc(data')
        ylim([0 length(data)])
        xlim([0 length(data)])
        ax = gca;
        ax.DataAspectRatio = [1 1 1];
        hold off
end

