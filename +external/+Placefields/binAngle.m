function binnedAngle=binAngle(at,gridsize)
%to be used in conjunction with binpos, so use the same gridsize.
%Change the time increment of time (second col) to gridsize
%interpolate angle (first col) for those time points.

tbin = linspace(at(1,2),at(end,2),gridsize);

% unwrap to prevent replication of the endpoint when interpolating polar
% coordinates, which can cause a sudden jump at 180 deg

theta_unwrapped = unwrap(at(:,1));
binnedAngle = interp1(at(:,2),theta_unwrapped,tbin); %still unwrapped

