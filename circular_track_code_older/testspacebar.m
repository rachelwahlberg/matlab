%testspacebar
% how to exit out of while loop using space bar


f = figure('Name','Keep me in focus so I can respond to space bar');
set(f,'KeyPressFcn',@isspacebar); %see below

drawnow;  %important
%get(f,'UserData')
i =0;
while i == 0
    drawnow; %important
end

disp('space bar pushed')
disp(['i = ' num2str(i)])

function isspacebar(src,event)
if strcmp(event.Key,'space')
    disp('happy')
    assignin("caller","i",1)
    i = 1;
end

end
