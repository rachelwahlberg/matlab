% Get fig one
open('14-Jun-2023_12.20.07.380_Harry_session6_performance.fig')
h = get(gcf,'Children'); % to get subplots

overallScatter = findobj(h(2),'Type','scatter');
overallXdata{1} = [overallScatter(:).XData];
overallYdata{1} = [overallScatter(:).YData];

ioScatter = findobj(h(1),'Type','Scatter');

i = 1;
j = 1;
for h = 1:length(ioScatter)
    if strcmp(ioScatter(h).MarkerFaceColor,'none') == true
w(i) = h;
i = i + 1;
    else
        b(j) = j;
        j = j+1;
    end
end

wXdata{1} = [ioScatter(w).XData];
wYdata{1} = [ioScatter(w).YData];
bXdata{1} = [ioScatter(b).XData];
bYdata{1} = [ioScatter(b).YData];

clear h overallScatter ioScatter

% Get fig two
open('14-Jun-2023_01.06.14.139_Harry_session6_performance.fig')
h = get(gcf,'Children'); % to get subplots

overallScatter = findobj(h(2),'Type','scatter');
overallXdata{2} = [overallScatter(:).XData];
overallYdata{2} = [overallScatter(:).YData];

ioScatter = findobj(h(1),'Type','Scatter');


i = 1;
j = 1;
for h = 1:length(ioScatter)
    if strcmp(ioScatter(h).MarkerFaceColor,'none') == true
w(i) = h;
i = i + 1;
    else
        b(j) = j;
        j = j+1;
    end
end

wXdata{2} = [ioScatter(w).XData];
wYdata{2} = [ioScatter(w).YData];
bXdata{2} = [ioScatter(b).XData];
bYdata{2} = [ioScatter(b).YData];

%% make combined figure

figure
hold all
grid_height = 2; grid_width = 2;
h1 = tiledlayout(grid_height,grid_width);
position = 1; h = 1; w = 1; ax1 = nexttile(position,[h,w]); %position
position = 2; h = 1; w = 1; ax2 = nexttile(position,[h,w]); %raster
position = 3; h = 1; w = 1; ax3 = nexttile(position,[h,w]); %
position = 4; h = 1; w = 1; ax4 = nexttile(position,[h,w]); % acg whole session

axes(ax1)
hold on
title('Overall Performance')
xlabel('Minutes from Port 1 start')
plot(overallXdata{1},overallYdata{1},'k','LineWidth',2)
hold off

axes(ax2)
hold on
xlabel('Minutes from Port 2 start')
plot(overallXdata{2},overallYdata{2},'k','LineWidth',2)
hold off

linkaxes([ax1 ax2],'y')

axes(ax3)
hold on
plot(wXdata{1},wYdata{1},'cyan','LineWidth',2)
plot(bXdata{1},bYdata{1},'black','LineWidth',2)
xlabel('Minutes from Port 1 start')
hold off

axes(ax4)
hold on
plot(wXdata{2},wYdata{2},'cyan','LineWidth',2)
plot(bXdata{2},bYdata{2},'black','LineWidth',2)
xlabel('Minutes from Port 2 start')
hold off

linkaxes([ax3 ax4],'y')




















