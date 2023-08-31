function allHVS = getHighVoltageSignals(wlfp,corticalCh,fs)

% adapted from Uktu code, "HighVoltageSignals" 
% data should be sent in as lfp.data
% corticalCh and hpcCh as vector of channels
% put in either premazelfp, mazelfp, or postmazelfp

%%%%%%%%%%%%%%%%%%%%%%% Set params %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

params.Fs = fs;

%%%%%%%%%%%%%%%%%%%%%%%% Cortical %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for ch = 1:numel(corticalCh)
    data = wlfp(:,ch);
   
    params.fpass = [5 10];
    [Stheta] = mtspecgramc(data,[2,1],params);
    meanStheta_chs(ch,:) = mean(abs(Stheta),2,'omitnan');

    params.fpass = [20 40]; % 
    [Sabovetheta,t,f] = mtspecgramc(data,[2,1],params);
    meanSabovetheta_chs(ch,:) = mean(abs(Sabovetheta),2,'omitnan');
end

meanStheta = mean(meanStheta_chs,1,'omitnan');
meanSabovetheta = mean(meanSabovetheta_chs,1,'omitnan');

z_theta = (meanStheta-mean(meanStheta,'omitnan'))/std(meanStheta,'omitnan');
z_abovetheta = (meanSabovetheta-mean(meanSabovetheta,'omitnan'))/std(meanSabovetheta,'omitnan');

%as a check
figure
hold all
plot(t,z_theta,'r')
plot(t,z_abovetheta,'g')

idxTheta_cor = z_theta > 1.2;
idxAbovetheta_cor = z_abovetheta > 1.2;

clear data meanStheta_chs meanSabovetheta_chs

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% HPC %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% for ch = 1:numel(hpcCh) 
%     data = mazelfp(:,hpcCh(ch));
%    
%     params.fpass = [5 10];
%     [Stheta] = mtspecgramc(data,[2,1],params);
%     meanStheta_chs(ch,:) = mean(Stheta,1);
% 
%     params.fpass = [10 40]; % 
%     [Sabovetheta] = mtspecgramc(data,[2,1],params);
%     meanSabovetheta_chs(ch,:) = mean(Sabovetheta,1);
% end
% 
% clear data
% 
% z_theta = zscore(meanStheta_chs);
% z_abovetheta = zscore(meanSabovetheta_chs);
% 
% %as a check
% figure
% hold all
% plot(t,z_theta)
% plot(t,z_abovetheta)
% 
% idxTheta_hpc = z_theta > 2;
% idxAbovetheta_hpc = z_abovetheta > 2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Combine %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

HVSevents = zeros(size(z_theta));
HVSevents = double(idxAbovetheta_cor);
%HVSevents(HVSevents == 0) = nan;

st = 1;
sp = 1;
for i = 2:length(HVSevents)
    if HVSevents(i)-HVSevents(i-1) == 1
        firstpass_start(st) = t(i);
        st = st +1 ;
        continue
    end

    if HVSevents(i)-HVSevents(i-1) == -1
        firstpass_stop(sp) = t(i);
        sp = sp+1;
        continue
    end
end

firstpass(:,1) = firstpass_start;
firstpass(:,2) = firstpass_stop;

%%%%%%%%% combine if they're super close %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

minInterSignalPeriod = 5; % 1 second to start
hvs = firstpass(1,:);
secondpass = [];
for i = 2:size(firstpass,1)
	if firstpass(i,1) - hvs(2) < minInterSignalPeriod
		% Merge
		hvs = [hvs(1) firstpass(i,2)];
	else
		secondpass = [secondpass ; hvs];
		hvs = firstpass(i,:);
	end
end
secondpass = [secondpass;hvs];
allHVS = secondpass;

%%%%%%%%%%%%%%%%%%%%%% Plot check %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% figure
% hold all
% imagesc(t,f,sp) 

% hvs_start2 = hvs_start * 1250;
% hvs_stop2 = hvs_stop * 1250; % now in SAMPLES. 

% figure
% hold all
% plot(mazelfp(:,1))
% xline(hvs_start2,'g')
% xline(hvs_stop2,'r')
% 
% 















