function [abovethresh_nspks,goodindex,belowthresh_nspks,badindex] = velocityfilterKD(xyt,nspk,thresh)
%must be AT LEAST the threshold - eliminates slower values

%take the differential of x values and then y values to determine points
%for which speed is below the provided threshold

speed2d = sqrt()

badindex = find(abs(diff(xyt(:,1)))<thresh | abs(diff(xyt(:,2)))<thresh); %added in for 2D

badindex = find(abs(diff(xyt(:,1)))<thresh | abs(diff(xyt(:,2)))<thresh); %added in for 2D
goodindex = find(abs(diff(xyt(:,1)))>=thresh & abs(diff(xyt(:,2)))>=thresh);


belowthresh_nspks = zeros(size(nspk));
belowthresh_nspks(badindex) = nspk(badindex); %the below
nspk(badindex)=0;
nspk(badindex+1)=0; % just to be safe

abovethresh_nspks = nspk;

