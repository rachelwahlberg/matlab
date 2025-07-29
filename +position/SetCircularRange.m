function [xztnew, spkCountNew, raster, edges, spktimes,trial] = SetCircularRange(xzt,spkCount,limits)
%based off of KD setXRange - adapted for RW circular track.

% xyt is normalized x and y position (between 0 and 1), along with time. they are already binned
% into specific time bins.

% limits.x = [xOutsideLeft xInsideLeft xInsideRight xOutsideRight]
% limits.z = [zLowerBottom zLowerTop zUpperBottom zUpperTop]

% function [xytnew, nspknew, raster, edges, spktimes] = position.SetCircularRange(xyt,nspk,limits,doraster)
%
% function takes xyt and nspk, the spikes during xyt, and outputs only the
% xyt and nspk's within 'limits', which should correspond to indices of xyt
% also gives out for 'raster' the number of spikes versus x-position  
% along with the 'edges' defined by limits, the indices 

if (size(spkCount,1)>size(spkCount,2))
    spkCount = spkCount';
end

%if (nargin<4); doraster = 1; end

%preallocate
xztnew = xzt;spkCountNew = spkCount;raster = [];edges = [];spktimes = [];trial =[];

if isempty(limits)
     warning('Input limits of SetPositionRange are empty')
    return
end

limitsx = xzt(:,1) > limits.x(1) & xzt(:,1) < limits.x(2);
limitsz = xzt(:,2) > limits.z(1) & xzt(:,2) < limits.z(2);
innercircle = limitsx & limitsz;
innercircle = find(innercircle == 1);
xztnew(innercircle) = 0;
spkCountNew(innercircle) = 0;

% 
% 
% 
% 
% 
% 
% 
%     for ii = 1:length(limits(1,:))
%         xztnew = [xzt(limits(1,ii):limits(2,ii),:); xztnew];
%         nspknew = [nspk(limits(1,ii):limits(2,ii)) nspknew];
%         trial = [ii*ones(limits(2,ii)-limits(1,ii)+1,1); trial];
%         if doraster
%             spkind = find(nspk(limits(1,ii):limits(2,ii))) + limits(1,ii)-1 ;
%             raster = [raster; [ii*ones(length(spkind),1) xzt(spkind,1)]];
%             spktimes = [spktimes; xzt(spkind,size(xzt,2))];
%         end
%     end
% end
% 
% if doraster
%     edges = [1:length(limits(1,:));xzt(limits(1,:),1)';xzt(limits(2,:),1)'];   % these mark the range of x values that are included in raster
% end

% odd and even indicate edges
