changeindex = [];
for a = 1:length(sessions)
    for s = 1:length(sessions{a})
        fnum = sessions{a}(s);
        perf = animalStructures{a}.preTables{fnum}.Performance;

        if ~(a==2 && fnum ==17)
        [m,i] = min(perf); % taken to be the point that pair shifts. works except for minerva file 17
        changeindex(a,s) = i;

        else
            [~,i] = min(perf(50:end)); %cause switch is at second min here
        end


    end
end

