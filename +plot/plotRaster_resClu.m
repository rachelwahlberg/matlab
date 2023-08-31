function [] = plotRaster_resClu(data,timewindow)

% data in res in first column, clu in second column
% time window in s
data.res=data.res/30000; % into seconds

res = data.res;
clu = data.clu;

idx = res>timewindow(1) & res<timewindow(2);

newdata = data(idx,:);


uniqueClu = unique(newdata.clu);

% sort by FR
for i = 1:length(uniqueClu)
nspikes(i) =sum(newdata.clu == uniqueClu(i));
end
[val,idx] = sort(nspikes,'descend');

uniqueClu1=uniqueClu(idx);
hold on
for i= 1:length(uniqueClu1)
    c = uniqueClu1(i);
    cidx = newdata.clu == c;
    resVals = newdata.res(cidx);
    
    scatter(resVals,ones(length(resVals),1)*i,"|","k",'LineWidth',1.2);
end

xlim(timewindow)
xvals = round(linspace(0,timewindow(2)-timewindow(1),10));
xticklabels(xvals)
xlabel('time (s) from trial start')
ylabel('Cell number')

title('Spike Raster for trial 6')

hold off



end