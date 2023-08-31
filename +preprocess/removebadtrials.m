function trials = removebadtrials(rawtrials,behaviortimestamps)

% currently rawtrials is in seconds from start of LFP 

% now in ms from lfp clock start
removestart = lfptimestamps.fromlfpclockstart(removetimestamps.sampformat(:,1));
removeend = lfptimestamps.fromlfpclockstart(removetimestamps.sampformat(:,2));


% right now NOT removing pre or post periods that are within the bad
% timeperiods

keep = 1;
for tri = 1:length(rawtrials.startTrial) % 
    msstart = seconds(rawtrials.startTrial(tri))*1000;
    msend = seconds(rawtrials.endTrial(tri))*1000;

    for remove = 1:length(removestart)

        if msstart >
    if msstart > removestart(tri) && msstart < removeend(tri)
        continue
    end

    if msend  > removestart(tri) && msstart < removeend(tri)
        continue
    end

        end
    end
end
