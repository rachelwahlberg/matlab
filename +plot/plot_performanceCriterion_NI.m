function [] = plot_performanceCriterion_NI(thresh,starttime)
thresh = .85; % 85% performance
continoo = 1;
n = 0;
totalcorrect = 0;
starttime = [];

if isempty(starttime)
starttime = datetime;
end

figure
hold on
title(['Performance de ' rat])
xlabel('minutes from start')
ylabel('performance')
while continoo == 1
    prompt = "Correct? 1 for yes, 0 for no";
    outcome = logical(input(prompt)); 
    totalcorrect = totalcorrect + outcome; 
    n=n+1;
    performance = totalcorrect/n;
    d = datetime;
    scatter(minutes(d-starttime),performance)
    disp(['Performance = ' num2str(performance)])

    
    if performance >= 0.85 
       disp('YOU HAVE HIT THRESHOLD !! CELEBRATION !!')
       %saveas(gcf,[char(datetime("today")) '_' rat '_performance.fig'])
    end
end

% 1 = white, outside, correct


    