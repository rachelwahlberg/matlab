function [xytnew, nspknew, raster, edges, spktimes,trial] = SetXRange(xyt,nspk,limits,doraster)
% function [xytnew, nspknew, raster, edges, spktimes] = SetXRange(xyt,nspk,limits,doraster)
%
% function takes xyt and nspk, the spikes during xyt, and outputs only the
% xyt and nspk's within 'limits', which should correspond to indices of xyt
% also gives out for 'raster' the number of spikes versus x-position  
% along with the 'edges' defined by limits, the indices 

if (size(nspk,1)>size(nspk,2));
    nspk = nspk';
end

if (nargin<4); doraster = 1; end

xytnew = [];nspknew = [];raster = [];edges = [];spktimes = [];trial =[];
if isempty(limits)
%     warning('Input limits of SetXRange are empty')
    return
end

for ii = 1:length(limits(1,:));
        xytnew = [xyt(limits(1,ii):limits(2,ii),:); xytnew];
        nspknew = [nspk(limits(1,ii):limits(2,ii)) nspknew];
        trial = [ii*ones(limits(2,ii)-limits(1,ii)+1,1); trial];
        if doraster
            spkind = find(nspk(limits(1,ii):limits(2,ii))) + limits(1,ii)-1 ;
            raster = [raster; [ii*ones(length(spkind),1) xyt(spkind,1)]];
            spktimes = [spktimes; xyt(spkind,size(xyt,2))];
        end
end

if doraster
    edges = [1:length(limits(1,:));xyt(limits(1,:),1)';xyt(limits(2,:),1)'];   % these mark the range of x values that are included in raster
end

% odd and even indicate edges
