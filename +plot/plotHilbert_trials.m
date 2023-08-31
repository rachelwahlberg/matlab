function [] = plotHilbert_trials(lfp,trials,fs)

% for preliminary figure
% poitns 80 to 180 into trial (with 60 sr) =
rangestart = floor(80/60*1000);
rangeend = ceil(180/60*1000);

% narrow bandpass for theta and gamma
figure
hold on
xlabel('time from trial start (ms)')
xvals = linspace(rangestart/1250,rangeend/1250,10)*1000;
xvals = round(xvals,1);


for tri = 1:size(trials,1)
    title('Bandpassed theta and gamma for sample trial 5')
    tristart = trials.startTrial_file(tri);
    triend = trials.endTrial_file(tri);
    trialdata = lfp.data(tristart:triend,126); % representative from hpc ch

    wdata = WhitenSignal(trialdata);
    fs.narrowhi = 5;
    fs.narrowlo = 12;
    [~,theta] = htfilter(wdata,fs);

    fs.narrowhi = 50;
    fs.narrowlo = 70;
    [~,gamma] = htfilter(wdata,fs);

    
    plot(real(gamma(rangestart:rangeend)),'red','LineWidth',1);
    plot(real(theta(rangestart:rangeend)),'blue','LineWidth',1);

    xticklabels(xvals);
    xlim([1 1668])
    legend([{'50-70Hz'},{'5-12Hz'}])

    pause
    clf
    hold on
end

end
