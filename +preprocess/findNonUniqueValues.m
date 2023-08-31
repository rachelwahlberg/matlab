% finding incorrect repeat values in behaviortimestamps.fromlfpfilestarts
function idx = findNonUniqueValues(timestamps)

u=unique(timestamps.fromlfpfilestart);         % the unique values
[n,bin]=histc(timestamps.fromlfpfilestart,u);  % count how many of each and where
%ix1=find(n>1);       % index to bin w/ more than one

idx=[];
for v=find(n>1).'
  idx=[idx;find(bin==v)]; % idx goes back into the original
end

end

% 34220:34349

