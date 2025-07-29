function [cellIdx,niIdx] = convertSensorName(pcbIdx)

% to get PCB board sensor number into the convention used on the
% national instruments

% the commented out cellIdx is the original conversion if you put into
% digital input [0 1 2 3 4 5 6 ...]; 

for i = 1:length(pcbIdx)
    if pcbIdx(i) == 0
        cellIdx(i) = 1;
        niIdx{i} = 'port0/line5';
        %cellIdx(i) = 6;  
        continue
    elseif pcbIdx(i) == 1
        cellIdx(i) = 2;
        niIdx{i} = 'port0/line1';
        %cellIdx(i) = 2;
        continue
    elseif pcbIdx(i) == 2
        cellIdx(i) = 3;
        niIdx{i} = 'port0/line6';
        %cellIdx(i) = 7;
        continue
    elseif pcbIdx(i) == 3
        cellIdx(i) = 4;
        niIdx{i} = 'port0/line2';
        %cellIdx(i) = 3;
        continue
    elseif pcbIdx(i) == 4
        cellIdx(i) = 5;
        niIdx{i} = 'port0/line7';
        %cellIdx(i) = 8;
        continue
    elseif pcbIdx(i) == 5
        cellIdx(i) = 6;
        niIdx{i} = 'port0/line3';
        %cellIdx(i) = 4;
        continue
    elseif pcbIdx(i) == 6
        cellIdx(i) = 7;
        niIdx{i} = 'port2/line3';
        %cellIdx(i) = 20;
        continue
    elseif pcbIdx(i) == 7
        cellIdx(i) = 8;
        niIdx{i} = 'port2/line2';
        %cellIdx(i) = 19;
        continue
    elseif pcbIdx(i) == 8
        cellIdx(i) = 9;
        niIdx{i} = 'port2/line6';
        %cellIdx(i) = 23;
        continue
    elseif pcbIdx(i) == 9
        cellIdx(i) = 10;
        niIdx{i} = 'port2/line4';
        %cellIdx(i) = 21;
        continue
    elseif pcbIdx(i) == 10
        cellIdx(i) = 11;
        niIdx{i} = 'port2/line1';
        %cellIdx(i) = 18;
        continue
    elseif pcbIdx(i) == 11
        cellIdx(i) = 12;
        niIdx{i} = 'port2/line0';
        %cellIdx(i) = 17;
        continue
    elseif pcbIdx(i) == 12
        cellIdx(i) = 13;
        niIdx{i} = 'port1/line7';
        %cellIdx(i) = 16;
        continue
    elseif pcbIdx(i) == 13
        cellIdx(i) = 14;
        niIdx{i} = 'port1/line6';
        %cellIdx(i) = 15;
        continue
    elseif pcbIdx(i) == 14
        cellIdx(i) = 15;
        niIdx{i} = 'port2/line7';
        %cellIdx(i) = 24;
        continue
    elseif pcbIdx(i) == 15
        cellIdx(i) = 16;
        niIdx{i} = 'port1/line5';
        %cellIdx(i) = 14;
        continue
    end

end








