function length = getfilelength(kilobytes)

nch = 70;
samplingrate = 30000;
nbits = 16;
m=1;
length = kilobytes*1024/(samplingrate*nbits*nch/8); % gives length in seconds

if m == 1
    length = length/60; % to give output in min
end