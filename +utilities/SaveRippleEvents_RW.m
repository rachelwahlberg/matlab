function SaveRippleEvents_RW(filename,ripples,channelID)

%SaveRippleEvents - Save hippocampal ripple (~200Hz oscillations) events.
%
%  USAGE
%
%    SaveRippleEvents(filename,ripples,channelID)
%
%    filename       file to save to
%    ripples        ripple info as provided by <a href="matlab:help FindRipples">FindRipples</a>
%    channelID      channel ID (appended to the event description)
%
%  SEE
%
%    See also FindRipples, RippleStats, PlotRippleStats, SaveEvents.

% Copyright (C) 2004-2011 by Michaël Zugaro
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

if nargin < 3,
  error('Incorrect number of parameters (type ''help <a href="matlab:help SaveRippleEvents">SaveRippleEvents</a>'' for details).');
end
% 
% n = size(ripples,1);
% r = ripples(:,1:3)';
% events.time = r(:);
% for i = 1:3:3*n,
% 	events.description{i,1} = ['Ripple start ' int2str(channelID)];
% 	events.description{i+1,1} = ['Ripple peak ' int2str(channelID)];
% 	events.description{i+2,1} = ['Ripple stop ' int2str(channelID)];
% end
% 
% SaveEvents_RW(filename,events);

%%%%%%

%n = size(ripples,1);
%r = ripples(:,1:3)';
% events.time = r(:);
% for i = 1:3:3*n,
% 	events.description{i,1} = ['Ripple start ' int2str(channelID)];
% 	events.description{i+1,1} = ['Ripple peak ' int2str(channelID)];
% 	events.description{i+2,1} = ['Ripple stop ' int2str(channelID)];
% end

r = [ripples.timestamps(:,1),ripples.peaks(:,1),ripples.timestamps(:,2)]';
r = r*1000;
n = size(r,2)*3;
% %working version
% for i = 1:3:n
%     events{i,1} = r(i);
%     events{i+1,1} = r(i+1);
%     events{i+2,1} = r(i+2);
% 
%     events{i,2} = 'ripple start';
%     events{i+1,2} = 'ripple peak';
%     events{i+2,2} = 'ripple stop';
% end
%but getting into format for SaveEvents func
for i = 1:3:n
    events.time(i,1) = r(i);
    events.time(i+1,1) = r(i+1);
    events.time(i+2,1) = r(i+2);

    events.description{i} = ['ripple start chan ',int2str(channelID)];
    events.description{i+1} = ['ripple peak chan ',int2str(channelID)];
    events.description{i+2} = ['ripple stop chan ',int2str(channelID)];
end

SaveEvents_RW(filename,events);