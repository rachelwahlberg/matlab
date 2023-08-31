%% Get zIdphi (z score integrated change in Phi) for each trial
% using Papale's code (odorTrails on his github)

%% Load position data
%can make more versatile in the future

% get timestamps relative to first lfp file
load(['/home/wahlberg/Exp_Data/M1_Nov2021/20211123/merged_M1_20211123_raw/' ...
    'merged_M1_20211123_raw_60Hz_behaviortimestamps.mat']); %behavior_timestamps
% get optitrack position data
load(['/home/wahlberg/Exp_Data/M1_Nov2021/20211123/merged_M1_20211123_raw/' ...
    'merged_M1_20211123_raw.dat.optitrack.behavior.mat']); %optitrack
% get starttime of first lfp file
load(['/home/wahlberg/Exp_Data/M1_Nov2021/20211123/merged_M1_20211123_raw/' ...
    'merged_M1_20211123_raw_1250Hz_lfptimestamps.mat']); %lfp_timestamps

% check that it seems accurate so far; each trial a new color
figure
hold all
for t= 1:ntrials
    plot(x{t},y{t})
end
xdecision = -10000;
ydecision = [-3000 3000]; 
trialdata = getYmazeTrials()
%% for each trial - 
% Also think about whether the length of the trial causes any difference?
% Make sure you totally understand how the value itself is scored

allidphi = []; % will become the integration over the entire 
for t = 1:ntrials
    if t == 4
        idphi{t} = nan;
        idphi_mean(t) = nan;
        continue
    end
    xdata = (trialdata.xdata{t})';
    ydata = (trialdata.ydata{t})';

    [idphi0,xinterp,yinterp] = AnalyzeOdorTrail_RW(xdata,ydata);

    trialdata.xinterp{t}=xinterp; %interpolated version of x
    trialdata.yinterp{t}=yinterp;
    allidphi = cat(1,allidphi,idphi0);
    idphi{t}=idphi0;
    idphi_mean(t) = mean(idphi{t});
end

n = length(find(~isnan(idphi_mean)));
z = (idphi_mean - mean(idphi_mean,'omitnan'))./(std(idphi_mean,'omitnan')/sqrt(n));

z = (x-mu/(sd/sqrt(n)));


zidphi_means=nanzscore(idphi_mean);
zidphi=nanzscore(abs(allidphi));

    maxz(t)=max(zIdPhi{t});
    meanz(t) = nanmean(zIdPhi{t});
figure
hold all
for t = 1:ntrials
    plot(zIdPhi{t})
end

figure
scatter(1:length(maxz),maxz)


y = (x - nanmean(x))./nanstd(x);

%%


figure
hold all
for t = 1:ntrials
    plot(trialdata.xinterp{t},trialdata.yinterp{t})
end

















