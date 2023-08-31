function [] = renamefiles(basepath)
%basepath = '/home/wahlberg/Analysis/M1_Nov2021/20211123/merged_M1_20211123/';
d= dir(basepath);
for f = 1:length(d)
    fname = d(f).name;
    index = strfind(fname,'raw');
    if index
        newname = [fname(1:index-2),fname(index+3:end)];
        movefile([basepath fname],[basepath newname])
    end
end
