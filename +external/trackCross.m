
function [outgoing,incoming] = trackCross(xx, xmin, xmax,movement);

%  function [outgoing, incoming] = trackCross(xx, xmin, xmax);
% this function returns the indices of xx for which the trajectory is outgoing,
% between xmin and xmax, or incoming, between xmax and xmin, while
% eliminating instances when the rat turned around halfway through.

if (nargin<4); movement = 1; end  % this is a flag to check if the animals movement is constant
    % It's optional, in case theta criterion is not used.

cross1 = xx >= xmin ;
cross2 = xx > xmax ;
      
cr1plus = find(diff(cross1)==+1);
cr1min =  find(diff(cross1)==-1);

cr2plus = find(diff(cross2)==+1);
cr2min =  find(diff(cross2)==-1);

outgoing = []; incoming = [];
badcrs = 0; % to be used later
for kk = 1:length(cr2plus);  %% for every instance in which 2nd gate is crossed
    ii = find((cr1plus < cr2plus(kk)),1,'last');  % find previous instance of 1st gate crossed
    if cr1plus(ii)~=badcrs;  % record
        if movement
            stoppoints = [cr1plus(ii); find(diff(xx(cr1plus(ii):cr2plus(kk))) < 0) + cr1plus(ii) ; cr2plus(kk)]; % 
            ilongest = find(diff(stoppoints)==max(diff(stoppoints)));
                    % find longest continuous movement along trajectory
            outgoing = [outgoing [stoppoints(ilongest)+1 ; stoppoints(ilongest+1)]]; % track validity times
        else
            outgoing = [outgoing [cr1plus(ii)+1 ; cr2plus(kk)]];
        end
        badcrs = cr1plus(ii);  % make sure this cross time is no longer used
    end
end
badcrs = 0;
for kk = 1:length(cr1min);  %% for every instance in which 1st gate is crossed
    ii = find((cr2min < cr1min(kk)),1,'last');  % find previous instance of 2nd gate crossed
    if cr2min(ii)~=badcrs;  
        if movement
            stoppoints = [cr2min(ii); find(diff(xx(cr2min(ii):cr1min(kk))) > 0) + cr2min(ii); cr1min(kk)]; % 
            ilongest = find(diff(stoppoints)==max(diff(stoppoints)));
            incoming = [incoming [stoppoints(ilongest)+1 ; stoppoints(ilongest+1)]]; % track validity times      
                    % make sure animal keeps 
%         cr2min(ii) = max([find(diff(xx(cr2min(ii):cr1min(kk))) > minimumdiff)+cr2min(ii);cr2min(ii)]);
%         incoming = [incoming [cr2min(ii);cr1min(kk)]]; % track validity
%         times
        else
            incoming = [incoming [cr2min(ii)+1 ; cr1min(kk)]]; % track validity
        end
        badcrs = cr2min(ii);  % make sure this cross time is no longer used
    end
end
