% PFPlot(PlaceMap, OccupancyMap, TopRate, TimeThresh, Gamma) =
%
% Plots a place field (see PFClassic)
%
% These are now plotted with directions as in the wheel file!
%
% NB PlaceMap and OccupancyMap should be in Hz and seconds.
% this is NOT how they are output from PFClassic (which depends)
% on the sample rates.  So convert them!


function PFPlot_RW(PlaceMap, sTimeSpent, TopRate, TimeThresh, Gamma)

FireRate = PlaceMap;
sTimeSpent = sTimeSpent;
FireRate = FireRate./max(FireRate,[],2);

if nargin<4
    TimeThresh=.1;
end
if nargin<5;
    Gamma = 1;
end
if nargin<3 | isempty(TopRate)
    TopRate = round(min([max(FireRate(find(sTimeSpent>TimeThresh))) 140]));
    if TopRate<1 , TopRate=1; end;
    if isempty(TopRate), TopRate=1; end;
end

Hsv(:,:,1) = (2/3) - (2/3)*external.clip(FireRate'/TopRate,0,1).^Gamma; %hue component determined by normalized FR w/ gamma correction
%2/3 stuff is to make it blue
Hsv(:,:,3) = 1./(1+TimeThresh./(sTimeSpent'+eps)); %value component is calculated based on time spent
Hsv(:,:,2) = ones(size(FireRate')); %saturation component

rgb = hsv2rgb(Hsv);
%nanidx = rgb(:,:,1) == 0;
%rgb(nanidx,:) = nan;
image(rgb); %function converts hue-saturation-value colors to rgb.
%set(gca, 'ydir', 'normal')

xlabels = round(linspace(0,2*pi,11),2);
% % most annoying bit is colorbar
h = gca;
h.YTickLabel = xlabels(2:end);
%h.YLabel = 'Track Position (radians)';
%h.DataAspectRatio=[1 1 1];

h2 = external.SideBar;
BarHsv(:,:,1) = (2/3) - (2/3)*(0:.01:1).^Gamma;
BarHsv(:,:,2) = ones(101,1);
BarHsv(:,:,3) = ones(101,1);
bar=(0:.01:1)*TopRate;
image(bar,0, hsv2rgb(BarHsv));
set(gca, 'ydir', 'normal');
set(gca, 'yticklabel',[]);
set(gca, 'xaxislocation', 'bottom');
xt = get(gca,'XTick');
if TopRate>max(xt)
    set(gca,'XTick',[xt TopRate])
end
axes(h);
