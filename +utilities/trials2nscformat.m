function events = trials2nscformat(trials)

% to convert RW trials format into format to put into neuroscope

% get into events.time and events.description

tri(1,:) = (trials.startTrial_file)'; % into double
tri(2,:) = (trials.endTrial_file)';

events.time = single(vertcat(tri(:))); % switches off trialstart, trialend

events.description = cell(length(events.time),1);
for i = 1:2:length(events.description)
    events.description{i} = "TrialStart";
    events.description{i+1} = "TrialStop";
end








end