function niIdx = convertSensorName(pcbIdx)

% to get PCB board sensor number into the convention used on the
% national instruments

for i = 1:length(pcbIdx)
    if pcbIdx(i) == 0
        niIdx{i} = 'port0/line5';
        continue
    elseif pcbIdx(i) == 1
        niIdx{i} = 'port0/line1';
        continue
    elseif pcbIdx(i) == 2
        niIdx{i} = 'port0/line6';
        continue
    elseif pcbIdx(i) == 3
        niIdx{i} = 'port0/line2';
        continue
    elseif pcbIdx(i) == 4
        niIdx{i} = 'port0/line7';
        continue
    elseif pcbIdx(i) == 5
        niIdx{i} = 'port0/line3';
        continue
    elseif pcbIdx(i) == 6
        niIdx{i} = 'port2/line3';
        continue
    elseif pcbIdx(i) == 7
        niIdx{i} = 'port2/line2';
        continue
    elseif pcbIdx(i) == 8
        niIdx{i} = 'port2/line6';
        continue
    elseif pcbIdx(i) == 9
        niIdx{i} = 'port2/line4';
        continue
    elseif pcbIdx(i) == 10
        niIdx{i} = 'port2/line1';
        continue
    elseif pcbIdx(i) == 11
        niIdx{i} = 'port2/line0';
        continue
    elseif pcbIdx(i) == 12
        niIdx{i} = 'port1/line7';
        continue
    elseif pcbIdx(i) == 13
        niIdx{i} = 'port1/line6';
        continue
    elseif pcbIdx(i) == 14
        niIdx{i} = 'port2/line7';
        continue
    elseif pcbIdx(i) == 15
        niIdx{i} = 'port1/line5';
        continue
    end

end








