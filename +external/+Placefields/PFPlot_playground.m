% PFPlot(PlaceMap, OccupancyMap, TopRate, TimeThresh, Gamma) =
%
% Plots a place field (see PFClassic)
%
% These are now plotted with directions as in the wheel file!
%
% NB PlaceMap and OccupancyMap should be in Hz and seconds.
% this is NOT how they are output from PFClassic (which depends)
% on the sample rates.  So convert them!


function PFPlot_playground(PlaceMap, sTimeSpent, TopRate, TimeThresh, Gamma)

FireRate = PlaceMap;

if nargin<4
    TimeThresh=0.1;
end
if nargin<5;
    Gamma = 1;
end
if nargin<3 | isempty(TopRate)
    TopRate = round(min([max(FireRate(find(sTimeSpent>TimeThresh))) 140]));
    if TopRate<1 , TopRate=1; end;
    if isempty(TopRate), TopRate=1; end;
end

Hsv(:,:,1) = 2/3 - (2/3)*external.clip(FireRate'/TopRate,0,1).^Gamma; %hue component determined by normalized FR w/ gamma correction
%2/3 stuff is to make it blue
Hsv(:,:,2) = ones(size(FireRate')); %saturation component
Hsv(:,:,3) = 1./(1+TimeThresh./(sTimeSpent'+eps)); %value component is calculated based on time spent

rgb = hsv2rgb(Hsv);
%nanidx = rgb(:,:,1) == 0;
%rgb(nanidx,:) = nan;
%image(rgb,AlphaDataMapping='scaled',AlphaData=alpha1); %function converts hue-saturation-value colors to rgb.
%set(gca, 'ydir', 'normal')

%%%%%%%%%%%%%%%%% TRYING TO INTRODUCE VIRIDIS COLORMAP %%%%%%% 
% % Map the RGB values to the colormap 
% % not working cause only interpolating blue dimension of rgb
% cmap = external.viridis_whiteNans;
% indexedImage = round(interp1(linspace(0, 1, size(cmap, 1)), ...
%     1:size(cmap, 1), rgb(:,:,3), 'linear', 'extrap')); 
% 
% % Apply the colormap 
% 
% coloredImage2 = ind2rgb(indexedImage2,cmap2); 
% 
% 
% % Display the result 
% 
% imshow(coloredImage); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% % most annoying bit is colorbar
h = gca;
h.DataAspectRatio=[1 1 1];

h2 = external.SideBar;
BarHsv2(:,:,1) = ((2/3) - (2/3)*(0:.01:1).^Gamma);
BarHsv2(:,:,2) = ones(101,1);
BarHsv2(:,:,3) = ones(101,1);
bar=(0:.01:1)*TopRate;
image(bar,0, hsv2rgb(BarHsv));
set(gca, 'ydir', 'normal');
set(gca, 'ytick', []);
set(gca, 'xaxislocation', 'bottom');
xt = get(gca,'XTick');
if TopRate>max(xt)
    set(gca,'XTick',[xt TopRate])
end
axes(h);
