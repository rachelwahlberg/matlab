% performance criterion for choosing inside or outside
% on black versus white side
rat = 'Harry';

thresh = .85; % 85% performance
continoo = 1;
bn = 0;
wn = 0;

b_totalcorrect = 0;
w_totalcorrect = 0;
starttime = [];

if isempty(starttime)
starttime = datetime;
end

figure
hold on
title(['Inside/Outside for ' rat])
xlabel('minutes from start')
ylabel('performance')
while continoo == 1
    prompt = "1 for b in, 2 for b out, 3 for w out, 4 for w in";
    outcome = input(prompt); 

    if outcome == 1 % 
        bn = bn + 1;
        b_totalcorrect = b_totalcorrect + 1;
        b_performance = b_totalcorrect/bn;
    elseif outcome == 2
        bn = bn + 1;
        b_performance = b_totalcorrect/bn;
    elseif outcome == 4
        wn = wn + 1;
        w_totalcorrect = w_totalcorrect + 1;
        w_performance = w_totalcorrect/wn;
    elseif outcome == 3
        wn = wn + 1;
        w_performance = w_totalcorrect/wn;
    end

    d = datetime;

    if outcome == 1 || outcome == 2
        scatter(minutes(d-starttime),b_performance,'black','filled')
        disp(['Black Performance = ' num2str(b_performance)])
    elseif outcome == 3 || outcome == 4
        scatter(minutes(d-starttime),w_performance,'black')
        disp(['White Performance = ' num2str(w_performance)])
    end

end

% 1 = white, outside, correct