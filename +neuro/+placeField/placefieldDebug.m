%test case,from within MethodsFigure.plotPlaceCells1 AND plotPlaceCellsKD
%unit =6;
class = 1;
smooth = 1;
speedThresh = 10; %cm/s

figure
grid_height = 2; grid_width = 8;

for unit = 1:length(su)

clf
h1 = tiledlayout(grid_height,grid_width);

%%% pulling in from KD code to allow changing units

spks = double(su(unit).TimesInSamples.SpikeTimes)/spikeSR; % get into seconds
spkCount2 = histcounts(spks,histedges)';
spkCount2(end)=0; spkCount2(1)=0; % find number of spikes per unit time.
%spkCount2 = external.velocityfilterUK_2D(binnedPos,spkCount2,10);

%%%%%%%%

title(h1,['Class ' num2str(class) ', unit ' num2str(unit)])
position = 1; h = 2; w = 1; ax1 = nexttile(position,[h,w]);
position = 2; h = 1; w = 1; ax2 = nexttile(position,[h,w]); 
position = 3; h = 1; w = 1; ax3 = nexttile(position,[h,w]);
position = 4; h = 1; w = 1; ax4 = nexttile(position,[h,w]); 
position = 5; h = 1; w = 1; ax5 = nexttile(position,[h,w]); 
position = 6; h = 1; w = 1; ax6 = nexttile(position,[h,w]); 
position = 7; h = 1; w = 1; ax7 = nexttile(position,[h,w]); 
position = 8; h = 1; w = 1; ax8 = nexttile(position,[h,w]); 
position = 10; h = 1; w = 1; ax10 = nexttile(position,[h,w]);
position = 11; h = 1; w = 1; ax11 = nexttile(position,[h,w]); 
position = 12; h = 1; w = 1; ax12 = nexttile(position,[h,w]); 
position = 13; h = 1; w = 1; ax13 = nexttile(position,[h,w]);
position = 14; h = 1; w = 1; ax14 = nexttile(position,[h,w]);
position = 15; h = 1; w = 1; ax15 = nexttile(position,[h,w]);
position = 16; h = 1; w = 1; ax16 = nexttile(position,[h,w]);

% AX1: plots position data vertically over time, and then the locations of unit firing
axes(ax1)
%title('Position over Time')
try
    sut{1}{unit}.plotOnTrack3D;
    %hold on
    %sut{c}{unit}.plotPlaneonTrack3D(obj.Events.ruleSwitch)
    %scatter(-50,obj.Event.ruleSwitch,0)
catch 
end

%%%% AJ 2D position %%%%%
axes(ax2)
hold on
title(' AJ 2D position')
plot(frm{1}{unit}.SpikeUnitTracked.PositionData.data.X,...
    frm{1}{unit}.SpikeUnitTracked.PositionData.data.Z,'b')
plot(sut{1}{unit}.TimesInSamples.X,sut{1}{unit}.TimesInSamples.Z,'.r')
ax2.DataAspectRatio = [1 1 1];
hold off

%%%% KD 2D position %%%%
axes(ax10)
hold on
title('KD 2D position')
plot(x,z,'y')
plot(sut{1}{unit}.TimesInSamples.X,sut{1}{unit}.TimesInSamples.Z,'.r');
ax10.DataAspectRatio = [1 1 1];
hold off

if smooth == 1

    %%%% AJ Occupancy %%%%%
    axes(ax3)
    hold on
    title('imagesc(frm.OccupancyMapSmooth)')
    imagesc(frm{1}{unit}.OccupancyMapSmooth')
    xlim([1 length(frm{1}{unit}.OccupancyMapSmooth)]);
    ylim([1 length(frm{1}{unit}.OccupancyMapSmooth)]);
    %ax2.DataAspectRatio = [1 1 1];
    hold off

    %%%% AJ Spike map %%%%%
    axes(ax4)
    hold on
    title('frm.SmoothedSpike')
    imagesc(frm{1}{unit}.SmoothedSpike')
    xlim([1 length(frm{1}{unit}.SmoothedSpike)]);
    ylim([1 length(frm{1}{unit}.SmoothedSpike)]);
    ax4.DataAspectRatio = [1 1 1];
    hold off

    %%%% AJ firing rate map imagesc %%%%%
    axes(ax5)
    hold on
    title('imagesc(frm.MapSmooth)')
    imagesc(frm{1}{unit}.MapSmooth')
    xlim([1 length(frm{1}{unit}.MapSmooth)]);
    ylim([1 length(frm{1}{unit}.MapSmooth)]);
ax5.DataAspectRatio = [1 1 1];
    hold off

    %%%% AJ firing rate map function %%%%%
    axes(ax6)
    hold on
    title('frm.plotSmooth')
    frm{1}{unit}.plotSmooth;
    xlim([1 length(frm{1}{unit}.OccupancyMapSmooth)]);
    ylim([1 length(frm{1}{unit}.OccupancyMapSmooth)]);
    ax6.DataAspectRatio = [1 1 1];
    hold off

    %%%% AJ place map imagesc %%%%%
    axes(ax7)
    hold on
    title('imagesc(pfm.MapSmooth)')
    imagesc(pfm{1}{unit}.MapSmooth')
    xlim([1 length(pfm{1}{unit}.MapSmooth)]);
    ylim([1 length(pfm{1}{unit}.MapSmooth)]);
    ax7.DataAspectRatio = [1 1 1];
    hold off

    %%%% AJ place map function %%%%%
    axes(ax8)
    hold on
    title('pfm.PlotSmooth')
    pfm{1}{unit}.plotSmooth;
    xlim([1 100]);
    ylim([1 100]);
    ax8.DataAspectRatio = [1 1 1];
    hold off

    %%%% KD place map %%%%%
    axes(ax16)
    hold on
    title('KD place map')
    [pc3,TimeSpent,nSpikes,sTimeSpent,snSpikes,pc1,pc2] = ...
    external.Placefields.PFClassic(binnedPos(:,1:2),spkCount2,smoothVal,nGrid,timebinsize);
    hold off

    %%%% KD Occupancy map %%%%%
    axes(ax11)
    hold on
    title('imagesc(sTimeSpent)')   
    imagesc(sTimeSpent')
    ylim([0 length(sTimeSpent)])
    xlim([0 length(sTimeSpent)])
    ax11.DataAspectRatio = [1 1 1];
    hold off

    %%%% KD Spike map %%%%%
    axes(ax12)
    hold on
    title('imagesc(sSpikes)')
    imagesc(snSpikes')
    ylim([0 length(snSpikes)])
    xlim([0 length(snSpikes)])
    ax12.DataAspectRatio = [1 1 1];
    hold off

    %%%% KD firing rate map %%%%%
    axes(ax13)
    hold on
    title('KD fire rate map')
    imagesc(pc1') %snSpikes./(sTimeSpent + eps)
    ylim([0 length(pc1)])
    xlim([0 length(pc1)])    
    ax13.DataAspectRatio = [1 1 1];

    hold off
    
    %%%% KD firing rate map %%%%%
    axes(ax14)
    hold on
    title('fire rate map, AJ filter')
    imagesc(pc2') %AJ filter,spk&occ smoothe
    ylim([0 length(pc2)])
    xlim([0 length(pc2)])
    ax14.DataAspectRatio = [1 1 1];
    hold off

    %%%% firing rate map, NK version %%%%%
    axes(ax15)
    hold on
    title('smooth (rawspk/smoothocc)')
    imagesc(pc3') %nat version
    ylim([0 length(pc3)])
    xlim([0 length(pc3)])
    ax15.DataAspectRatio = [1 1 1];
    hold off

else %plot raw

    %%%% AJ RAW Occupancy %%%%%
    axes(ax3)
    hold on
    title('imagesc(frm.OccupancyMap)')
    imagesc(frm{1}{unit}.OccupancyMap')
    xlim([1 length(frm{1}{unit}.OccupancyMap)]);
    ylim([1 length(frm{1}{unit}.OccupancyMap)]);
    ax3.DataAspectRatio = [1 1 1];
    hold off

    %%%% AJ RAW Spike map %%%%%
    axes(ax4)
    hold on
    title('Raw AJ Spike map')
    xlim([1 length(frm{1}{unit}.MapOriginal)]);
    ylim([1 length(frm{1}{unit}.MapOriginal)]);
    ax4.DataAspectRatio = [1 1 1];
    hold off

    %%%% AJ firing rate map %%%%%
    axes(ax5)
    hold on
    title('imagesc(frm.MapOriginal)')
    imagesc(frm{1}{unit}.MapOriginal')
    xlim([1 length(frm{1}{unit}.MapOriginal)]);
    ylim([1 length(frm{1}{unit}.MapOriginal)]);
    ax5.DataAspectRatio = [1 1 1];
    hold off

    %%%% AJ firing rate map function %%%%%
    axes(ax6)
    hold on
    title('frm.plot')
    frm{1}{unit}.plot;
      xlim([1 length(frm{1}{unit}.MapSmooth)]);
    ylim([1 length(frm{1}{unit}.MapSmooth)]);
    ax6.DataAspectRatio = [1 1 1];
    hold off

    %%%% AJ place map %%%%%
    axes(ax7)
    hold on
    title('imagesc(pfm.MapOriginal)')
    imagesc(pfm{1}{unit}.MapOriginal')
    xlim([1 length(pfm{1}{unit}.MapOriginal)]);
    ylim([1 length(pfm{1}{unit}.MapOriginal)]);
    ax7.DataAspectRatio = [1 1 1];
    hold off

    %%%% AJ place map function %%%%%
    axes(ax8)
    hold on
    title('pfm.plot')
    pfm{1}{unit}.plot;
    xlim([1 100]);
    ylim([1 100]);
    ax8.DataAspectRatio = [1 1 1];
    hold off

    %%%% KD place map %%%%%
    axes(ax14)
    hold on
    title('KD place map')
    [pc3,TimeSpent,nSpikes,sTimeSpent,snSpikes,pc1,pc2] = ...
    external.Placefields.PFClassic(binnedPos(:,1:2),spkCount2,smoothVal,nGrid,timebinsize);

    %%%% KD Occupancy map %%%%%
    axes(ax11)
    hold on
    title('imagesc(TimeSpent')
    imagesc(TimeSpent')
        ylim([0 length(TimeSpent)])
    xlim([0 length(TimeSpent)])
    ax11.DataAspectRatio = [1 1 1];
    hold off

    %%%% KD Spike map %%%%%
    axes(ax12)
    hold on
    title('imagesc(nSpikes)')
    imagesc(nSpikes')
        ylim([0 length(nSpikes)])
    xlim([0 length(nSpikes)])
    ax12.DataAspectRatio = [1 1 1];
    hold off

    %%%% KD firing rate map %%%%%
    axes(ax13)
    hold on
    title('smooth (rawspk/smoothocc)')
    imagesc(pc3') %nat version
    ylim([0 length(pc3)])
    xlim([0 length(pc3)])
    ax13.DataAspectRatio = [1 1 1];
    hold off
end

 %linkaxes([ax2 ax3 ax4 ax5 ax6 ax7 ax8 ax10 ax11 ax12 ax13 ax14 ax15 ax16],'xy')

 pause
 
end
