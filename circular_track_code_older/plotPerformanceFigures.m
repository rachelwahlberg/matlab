%plotPerformanceFigures

%add any new files for that rat to its folder%
foldername = 'E:\Rachel\Datapoints';
cd(foldername)

rat = 'Minerva';
movefile(['*' rat '*'], ['E:\Rachel\Datapoints\' rat])

%get filenames
animalDirectory = dir([foldername '\' rat '\*' rat '_pretraining*']);
filenames = {animalDirectory(:).name};

%formats for reading text files
lineFormat = '%s%s%d%f'; %string string integer double
columnNames = {'Date','Time','Port','Performance'};

%append to a file
for f = 1:length(filenames)
    data = readtable([foldername '\' rat '\' filenames{f}],'Format',lineFormat);
    data.Properties.VariableNames = columnNames;
   % perfTimes = data.Date;
   % perfTimes = datetime(strcat(data.Date,{' '},data.Time),"InputFormat","dd-mmm-yyyy hh:mm:ss:SSS");

    perfTimes = datetime(data.Time,"InputFormat","hh:mm:ss:SSS", "Format","HH:mm:ss:SSS");
    if any(perfTimes - perfTimes(1) < 0) %converting funny, it's making >noon be >midnight
       shiftthese = (perfTimes-perfTimes(1) < 0);
       perfTimes(shiftthese) = perfTimes(shiftthese) + hours(12);
    end %now it's correct time

    %days = table([],[],'VariableNames',{'relativeTimes','performance'})
    days(f).relativeTimes = minutes(perfTimes - perfTimes(1));
    days(f).performance = data.Performance;

  %  days = table('relativeTimes',relativeTimes,'performance', performance)
   
end

figure
hold on
xlabel('Time from zero (min)')
ylabel('Performance')
for f = 1:length(filenames)
plot(days(f).relativeTimes,days(f).performance)
end

figure
hold on
xlabel('Time from zero (min)')
ylabel('Performance')
plot(days(f).relativeTimes,days(f).performance,'LineWidth',2,'Color','k')
scatter(days(f).relativeTimes,days(f).performance,'Color','k')





