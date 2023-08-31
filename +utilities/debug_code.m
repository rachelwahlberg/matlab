
for c = 1:length(cells)
    for type = 1:length(spiketimes)
        for tri = 1:length(spiketimes{type})
            spk_index = find(spikeIDs{type}{tri} == cells(c));
            c_spiketimes{c}{type}{tri} =spiketimes{type}{tri}(spk_index);
        end
    end
end
