%% Eliminate Artifacts

% Get data
%dat = ft_fetch_data(data, 'header', hdr, 'begsample', trl(trlop,1)-fltpadding, 'endsample', trl(trlop,2)+fltpadding, 'chanindx', sgnind, 'checkboundary', strcmp(cfg.continuous, 'no'), 'skipcheckdata', 1);

lfp = bz_GetLFP('all','basepath',basepath); % or use m =memmapfile(merged_M1_20211123_raw.lfp','Format','int16') and then do the channels/division by SR manually
dat = lfp.data;

sumval = zeros(size(dat,1),1);
sum_sqr = zeros(size(dat,1),1);
numsmp = zeros(size(dat,1),1);

sumval = sumval + sum(dat,2,'omitnan'); 
sum_sqr = sum_sqr + sum(dat.^2,2,'omitnan');
numsmp = numsmp + size(dat,2); % number of samples

% compute the average and the standard deviation
datavg = sumval./numsmp;
datstd = sqrt(sum_sqr./numsmp - (sumval./numsmp).^2);

% 
smprate = 1250;
win_length =smprate*30; %30 seconds / window

numtrl=1; % FOR NOW, 6/29/22
zmax = cell(1, numtrl);
zsum = cell(1, numtrl);
zindx = cell(1, numtrl);

%preallocate z score matrices
zmax = zeros(size(dat));
zsum  = zeros(1,size(dat,1));
zindx = zeros(1,size(dat,1));

% get z scores    
    nsmp          = size(dat,1);
numwins = floor(nsmp/win_length);

for ch = 1:size(dat,2)
    for w = 1:numwins
        ind = 1;
        if w <= numwins
            win = ind:ind+win_length;
            zdata  = double(dat(win,ch)) - datavg(win)./datstd(win);  % convert the filtered data to z-values
            zsum   = sum(zdata,1,'omitnan');                  % accumulate the z-values over channels
            zmax(ch,w) = max(zdata,[],1);            % find the maximum z-value and remember it
        else
            win = ind:length(dat);
            zdata  = double(dat(win,ch)) - datavg(win)./datstd(win);  % convert the filtered data to z-values
            zsum   = sum(zdata,1,'omitnan');                  % accumulate the z-values over channels
            zmax(ch,w) = max(zdata,[],1);    
        end
    ind = w+win_length;   
    end
end

trlop = 1;
             zdata         = (dat - datavg(:,indvec(trlop)*ones(1,nsmp)))./datstd(:,indvec(trlop)*ones(1,nsmp));  % convert the filtered data to z-values
            zsum{trlop}   = nansum(zdata,1);                   % accumulate the z-values over channels
    [zmax{trlop},ind] = max(zdata,[],1);            % find the maximum z-value and remember it
    zindx{trlop}      = sgnind(ind); 
% create a vector that indexes the trials, or is all 1, in order
% to a per trial z-scoring, or use a static std and mean (used in lines 317
% and 328)
if pertrial
  indvec = 1:numtrl;
else
  indvec = ones(1,numtrl);
end
for trlop = 1:numtrl
  if strcmp(cfg.memory, 'low') % store nothing in memory (note that we need to preproc AGAIN... *yawn*
    fprintf('.');
    if hasdata
      dat = ft_fetch_data(data,        'header', hdr, 'begsample', trl(trlop,1)-fltpadding, 'endsample', trl(trlop,2)+fltpadding, 'chanindx', sgnind, 'checkboundary', strcmp(cfg.continuous, 'no'));
    else
      dat = ft_read_data(cfg.datafile, 'header', hdr, 'begsample', trl(trlop,1)-fltpadding, 'endsample', trl(trlop,2)+fltpadding, 'chanindx', sgnind, 'checkboundary', strcmp(cfg.continuous, 'no'), 'dataformat', cfg.dataformat);
    end
    dat = preproc(dat, cfg.artfctdef.zvalue.channel, offset2time(0, hdr.Fs, size(dat,2)), cfg.artfctdef.zvalue, fltpadding, fltpadding);
    zmax{trlop}  = -inf + zeros(1,size(dat,2));
    zsum{trlop}  = zeros(1,size(dat,2));
    zindx{trlop} = zeros(1,size(dat,2));
    
    nsmp          = size(dat,2);
    zdata         = (dat - datavg(:,indvec(trlop)*ones(1,nsmp)))./datstd(:,indvec(trlop)*ones(1,nsmp));  % convert the filtered data to z-values
    zsum{trlop}   = nansum(zdata,1);                   % accumulate the z-values over channels
    [zmax{trlop},ind] = max(zdata,[],1);            % find the maximum z-value and remember it
    zindx{trlop}      = sgnind(ind);                % also remember the channel number that has the largest z-value
  else
    % initialize some matrices
    zmax{trlop}  = -inf + zeros(1,size(dat{trlop},2));
    zsum{trlop}  = zeros(1,size(dat{trlop},2));
    zindx{trlop} = zeros(1,size(dat{trlop},2));
    
    nsmp          = size(dat{trlop},2);
    zdata         = (dat{trlop} - datavg(:,indvec(trlop)*ones(1,nsmp)))./datstd(:,indvec(trlop)*ones(1,nsmp));  % convert the filtered data to z-values
    zsum{trlop}   = nansum(zdata,1);                   % accumulate the z-values over channels
    [zmax{trlop},ind] = max(zdata,[],1);            % find the maximum z-value and remember it
    zindx{trlop}      = sgnind(ind);                % also remember the channel number that has the largest z-value
  end