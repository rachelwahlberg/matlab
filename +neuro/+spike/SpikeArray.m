classdef SpikeArray < neuro.spike.SpikeNeuroscope
    %SPIKEARRAY General spike analyses
    % not static; need to create a spike array object first.

    %%%%%%% Functions below: %%%%%%%%%%% 
%     obj = SpikeArray(spikeClusters,spikeTimes)
%     obj = getSpikeArrayWithAdjustedTimestamps(obj)
%     print(obj,varargin)
%     colstrings=tostring(obj,varargin)
%     obj=setTimeIntervalCombined(obj,ticd)
%     obj=setClusterInfo(obj,ci)
%     tbl=getspikeIDs(obj)
%     getMeanFireRateQuintiles(obj,nquintiles,timebininsec,order)
%     getMeanFireRate(obj,timebin)
%     []=plot(obj,tfm)
%     []=plotRaster(obj)
%     obj=getTimeInterval(obj,timeWindow)
%     spikeUnits=getSpikeUnits(obj,idx)
%     obj=get(obj,varargin)
%     pbe=getPopulationBurstEvents(obj) % EMPTY!
%     acg=getAutoCorrelogram(obj) % from neuro.spike
%     []=saveNeuroscopeFiles(obj,folder,filename)
%     [sa]=plus(obj,sa)
%     [obj]=setShank(obj,shankno)
%     [obj]=setLocation(obj,location)
%     [obj]=sort(obj,by)
%     st=getSpikeTimes(obj) %inherited
%     sc=getSpikeClusters(obj) % inherited

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    properties
        SpikeTable
        TimeIntervalCombined
        ClusterInfo
        Info
    end

    methods
        function obj = SpikeArray(spikeClusters,spikeTimes)
            %SPIKEARRAY Construct an instance of this class
            %   spiketimes should be in Timestamps.
            if ~isa(spikeClusters,'neuro.spike.SpikeArray')
                tablearray=horzcat( spikeTimes, double(spikeClusters));
                obj.SpikeTable=array2table(tablearray,'VariableNames',{'SpikeTimes','SpikeCluster'});
            else
                obj.SpikeTable=spikeClusters.SpikeTable;
                obj.TimeIntervalCombined=spikeClusters.TimeIntervalCombined;
                obj.ClusterInfo=spikeClusters.ClusterInfo;
                obj.Info=spikeClusters.Info;
            end
        end

        function obj = getSpikeArrayWithAdjustedTimestamps(obj)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            st=obj.SpikeTable;
            spiketimessample=st.SpikeTimes;
            ticd=obj.TimeIntervalCombined;
            adjustedspiketimessample=ticd.adjustTimestampsAsIfNotInterrupted(spiketimessample);
            obj.SpikeTable.SpikeTimes=adjustedspiketimessample;
        end
        function print(obj,varargin)
            logger=logging.Logger.getLogger;
            colstrings=obj.tostring(varargin{:});
            for icountgroup=1:numel(colstrings)
                colstring=colstrings{icountgroup};
                logger.info(colstring);
            end
        end
        function colstrings=tostring(obj,varargin)
            logger=logging.Logger.getLogger;
            tbl=obj.ClusterInfo;
            if nargin>1
                idx=ismember(varargin,tbl.Properties.VariableNames);
                if ~all(idx)
                    logger.warning([varargin{~idx} ' is not a column in the table.' strjoin(tbl.Properties.VariableNames,', ')]);
                end
                countGroup=varargin(idx);
            else
                countGroup={'group','ch'};
            end
            colstrings=cell.empty(numel(countGroup),0);
            for icountgroup=1:numel(countGroup)
                countgroup=countGroup{icountgroup};
                col=unique(tbl.(countgroup));
                colStr=cell.empty(numel(col),0);
                for icol=1:numel(col)
                    uniqueelement=col(icol);
                    idx=ismember(tbl.(countgroup),uniqueelement);
                    if iscategorical(uniqueelement)
                        uniqueelementstr=char(uniqueelement);
                    elseif isnumeric(uniqueelement)
                        uniqueelementstr=num2str(uniqueelement);
                    elseif isstring(uniqueelement)
                        uniqueelementstr=char(uniqueelement);
                    else
                        uniqueelementstr=uniqueelement;
                    end
                    colStr{icol}=[uniqueelementstr '(' num2str(sum(idx)) ')'];
                end
                colstring=strcat(countgroup,': ',strjoin(colStr,', '));
                colstrings{icountgroup}=colstring;
            end
        end
        function obj=setTimeIntervalCombined(obj,ticd)
            obj.TimeIntervalCombined=ticd;
        end
        function obj=setClusterInfo(obj,ci)
            obj.ClusterInfo=ci;
        end
        function tbl=getspikeIDs(obj)
            spikeIDs=unique(obj.SpikeTable.SpikeCluster);
            for iid=1:numel(spikeIDs)
                spikeid=spikeIDs(iid);
                spikecounts(iid)=sum(obj.SpikeTable.SpikeCluster==spikeid);
            end
            tbl=array2table(horzcat(spikeIDs,spikecounts'),'VariableNames',{'ID','count'});
        end
        %         function fr=getMeanFireRate(obj)
        %             sus=obj.getSpikeUnits;
        %             for isu=1:numel(sus)
        %                 su=sus(isu);
        %                 frs=su.getFireRate;
        %                 try
        %                     vals(isu,:)=frs.getValues;
        %                 catch
        %                 end
        %             end
        %             val_m=mean(vals,1);
        %             fr=Channel('Mean Over Units',val_m,frs.getTimeInterval);
        %         end
        function [frm, fre]=getMeanFireRateQuintiles(obj,nquintiles,timebininsec,order)
            sus=obj.getSpikeUnits;
            for isu=1:numel(sus)
                su=sus(isu);
                frs=su.getFireRate(timebininsec);
                try
                    vals1(isu,:)=frs.getValues;
                catch
                end
            end
            if ~exist('order','var')||isempty(order)
                vals=sort(vals1,1);
            else
                vals=vals1(order,:);
            end
            nunit=size(vals,1);
            qus=[round(quantile(1:nunit,nquintiles-1)) nunit];
            pre=1;
            for iquint=1:nquintiles
                idx=pre:qus(iquint);pre=qus(iquint)+1;
                thequint=vals(idx,:);
                themeanquint=mean(thequint,1);
                thesterrquint=std(thequint,1)/sqrt(size(thequint,1));
                frm{iquint}=neuro.basic.Channel(sprintf('Mean Over Units Quint, %d',iquint),...
                    themeanquint,frs.getTimeInterval);
                fre{iquint}=neuro.basic.Channel(sprintf('Mean Over Units Quint, %d',iquint),...
                    thesterrquint,frs.getTimeInterval);
            end
        end
        function [frm, fre]=getMeanFireRate(obj,timebin)
            sus=obj.getSpikeUnits;
            for isu=1:numel(sus)
                su=sus(isu);
                frs=su.getFireRate(timebin);
                try
                    vals(isu,:)=frs.getValues;
                catch
                end
            end
            themeanquint=mean(vals,1);
            thesterrquint=std(vals,1)/sqrt(size(vals,1));
            frm=neuro.basic.Channel(sprintf('Mean Over Units'),...
                themeanquint,frs.getTimeIntervalCombined);
            fre=neuro.basic.Channel(sprintf('Mean Over Units'),...
                thesterrquint,frs.getTimeIntervalCombined);

        end
        function frs=getFireRates(obj,timebininsec)
            sus=obj.getSpikeUnits;
            for isu=1:numel(sus)
                su=sus(isu);
                frs=su.getFireRate(timebininsec);
                try
                    vals(isu,:)=frs.getValues;
                catch
                end
            end
            frs=neuro.spike.FireRates(vals,[sus.Id],frs.getTimeIntervalCombined);
            frs.Info.TimebinInSec=timebininsec;
            frs.ClusterInfo=obj.ClusterInfo;
        end
        function []=plot(obj,tfm)

            angles1=angle(tfm.matrix);
            t=tfm.timeIntervalCombined.getTimePointsInAbsoluteTimes;
            spiketablet=obj.SpikeTable;
            ticd=obj.TimeIntervalCombined;
            clustinfo=obj.ClusterInfo;
            st=spiketablet.SpikeTimes;
            sc=spiketablet.SpikeCluster;
            figurename=sprintf('Rasterplot %s-%s',ticd.getRealTimeFor(st([1 end])));
            try close(figurename);catch,end
            figure('Name',figurename, 'Units','normalized','Position',[0 0 .2 .3])
            clustinfo1=clustinfo;
            locations=unique(clustinfo1.location);
            colors=linspecer(numel(locations),'qualitative');
            phasemap(8);
            color=phasemap(8);
            for iunit=1:height(clustinfo1)
                unit=clustinfo1(iunit,:);
                idx=ismember(sc,unit.id);
                stn=st(idx);
                if ~isempty(stn)
                    arr=ticd.getRealTimeFor(stn);
                    clear c
                    for ipoint=1:numel(arr)
                        pnt=arr(ipoint);
                        idx=find(t<pnt,1,'last');
                        angle1=angles1(idx);
                        idx_color=round((angle1+pi)/2/pi*(size(color,1)-1))+1;
                        c(ipoint,:)=color(idx_color,:);
                    end
                    hold on
                    idx_location=ismember(locations,unit.location);
                    s=scatter(arr,ones(size(arr))*iunit,ones(size(arr))*30,c,'|');
                    s.LineWidth=4;
                    s.MarkerEdgeAlpha=.8;
                    %                     p1=plot(arr,iunit...
                    %                         ,'Marker','|'...
                    %                         ,'LineWidth',2 ...
                    %                         ,'Color',colors(idx_location,:)...
                    %                         ,'MarkerSize',3,...
                    %                         'MarkerEdgeColor',colors(idx_location,:));
                    %                 s1=scatter(arr,ones(size(arr))*iunit,20,[0 0 0],'filled');
                    %                 s1.MarkerEdgeAlpha=.1;
                    %                 s1.MarkerFaceAlpha=.1;
                end
            end
            ax=gca;
            ax.YDir='reverse';
            ax.YTick=2:5:height(clustinfo1);
            ax.YTickLabel=clustinfo1.sh(ax.YTick);
            locs = find(diff(sign(angles1))>0);
            %             [~,locs]=findpeaks(angles1);
            t1=t(locs);
            for iline=1:numel(locs)
                xline(t1(iline));
            end
            phasebar('size',.15,'location','southeast');
        end
        function []=getRaster(obj,group,types) %group is which table row to sort by
            % GET raster here, plot it elsewhere (as subplot, or solo in plotRaster)

            spiketablet=obj.SpikeTable;

            clustinfo=obj.ClusterInfo;
            sortby = table2cell(clustinfo(:,group));
            sortby = string(sortby(:));

            if isempty(types) %hich types to use
            [types,~,gr]=unique(sortby); %'sh'
            end

            ticd=obj.TimeIntervalCombined;
            st=spiketablet.SpikeTimes;
            sc=spiketablet.SpikeCluster;
            colors=linspecer(numel(types),'qualitative');

            for ity=1:length(types)
                clustinfo1=clustinfo(gr == ity,:);  % just selecting for the units of that type

                t = string(ticd.getRealTimeFor(st([1 end])));
                figurename = strcat("Rasterplot_type",types(1), '_', t(1), "-", t(2));
                try close(figurename);catch,end
                figure('Name',figurename, 'Units','normalized','Position',[0 0 .2 .3])
                for iunit=1:height(clustinfo1)
                    unit=clustinfo1(iunit,:);
                    idx = ismember(sc,unit.id);
                    stn=st(idx);
                    if ~isempty(stn)
                        arr=ticd.getRealTimeFor(stn);
                        hold on
                       % idx_location=ismember(locations,unit.sh);
                        p1=plot(arr,iunit...
                            ,'Marker','|'...
                            ,'LineWidth',2 ...
                            ,'Color',colors(ity,:)...
                            ,'MarkerSize',3,...
                            'MarkerEdgeColor',colors(ity,:));
                        s1.MarkerEdgeAlpha=.1;
                        s1.MarkerFaceAlpha=.1;
                    end
                    iunit
                end
                ax=gca;
                ax.YDir='reverse';
                ax.YTick=2:5:height(clustinfo1);
                ax.YTickLabel=clustinfo1.sh(ax.YTick);
                ax.Title.String = figurename;

                if(figpath)
                saveas(ax,strcat(figpath,figurename, ".fig"));
                saveas(ax,strcat(figpath,figurename, ".pdf"));
                disp(strcat('Saved: ', figurename));
                end
                close(ax.Parent.Name);
            end
        end

        function [] = plotRaster(obj,group,types,figpath)



        end


        function obj=getTimeInterval(obj,timeWindow)
            if isdatetime(timeWindow)
                s=obj.TimeIntervalCombined.getSampleForClosest(timeWindow);
            elseif isduration(timeWindow)
                s=obj.TimeIntervalCombined.getSampleForClosest(obj.TimeIntervalCombined.getDate+timeWindow);
            end
            tbl=obj.SpikeTable;
            tbl((tbl.SpikeTimes<s(1))|(tbl.SpikeTimes>=s(2)),:)=[];
            tbl.SpikeTimes=tbl.SpikeTimes-s(1);
            obj.SpikeTable=tbl;
            obj.Info.TimeFrame=timeWindow;
            obj.TimeIntervalCombined=obj.TimeIntervalCombined.getTimeIntervalForTimes(timeWindow);
        end
        function spikeUnits=getSpikeUnits(obj,idx)
            tbl=obj.SpikeTable;
            ci=obj.ClusterInfo;
            if exist('idx','var') && ~isempty(idx)
                ci_sub=ci(idx,:);
            else
                ci_sub=ci;
            end
            for isid=1:height(ci_sub)
                aci=ci_sub(isid,:);
                spktimes=tbl.SpikeTimes(tbl.SpikeCluster==aci.id);
                spikeUnits(isid)=neuro.placecells.SpikeUnit(aci.id,spktimes,obj.TimeIntervalCombined);
                spikeUnits(isid)=spikeUnits(isid).setInfo(aci);
            end
        end
        function obj=get(obj,varargin)
            selected=true(height(obj.ClusterInfo),1);
            if nargin==1
            elseif nargin==2&& (isnumeric(varargin{1})||islogical(varargin{1}))
                selected=varargin{1};
            else
                cluinf=obj.ClusterInfo;
                stringinterest={'location','group'};
                selected_arg=false(height(obj.ClusterInfo),1);
                for iarg=1:nargin-1
                    arg=varargin{iarg};
                    for iint=1:numel(stringinterest)
                        interest=stringinterest{iint};
                        if any(ismember(cluinf.(interest),arg))
                            selected_arg=selected_arg|ismember(cluinf.(interest),arg);
                        end
                    end
                end
                selected=selected&selected_arg;
            end
            tbl=obj.SpikeTable;
            obj.ClusterInfo=obj.ClusterInfo(selected,:);
            obj.SpikeTable=tbl(ismember(tbl.SpikeCluster,obj.ClusterInfo.id),:);
        end
        function pbe=getPopulationBurstEvents(obj)

        end
        function acg=getAutoCorrelogram(obj)
            sus=obj.getSpikeUnits;
            acg=neuro.spike.AutoCorrelogram(sus);
        end
        function []=saveNeuroscopeFiles(obj,folder,filename)
            if ~exist('filename','var')
                filename='neuro';
            end
            if ~exist('folder','var')
                folder='.';
            end
            obj.saveCluFile(fullfile(folder,[filename  '.clu.0']));
            obj.saveResFile(fullfile(folder,[filename  '.res.0']));
        end
        function [sa]=plus(obj,sa)
            shift=max(obj.ClusterInfo.id);
            sa.ClusterInfo.id=sa.ClusterInfo.id+shift;
            sa.SpikeTable.SpikeCluster=sa.SpikeTable.SpikeCluster+shift;
            sa.ClusterInfo=sortrows([obj.ClusterInfo; sa.ClusterInfo],{'group','sh','ch'});

            sa.SpikeTable=sortrows([obj.SpikeTable; sa.SpikeTable],{'SpikeTimes'});
        end
        function [obj]=setShank(obj,shankno)
            obj.ClusterInfo.sh=ones(height(obj.ClusterInfo),1)*shankno;
        end
        function [obj]=setLocation(obj,location)
            str= repmat(location,height(obj.ClusterInfo),1);
            obj.ClusterInfo.location=str;
        end
        function [obj]=sort(obj,by)
            %by={'group','sh','ch'}
            obj.ClusterInfo=sortrows(obj.ClusterInfo,by);
        end

        function st = getUnitSpikeTimes(obj,id) % pull st for one unit
            st = table('SpikeTimes');
            idx = find(obj.SpikeTable.SpikeCluster == id);
            st = obj.SpikeTable.SpikeTimes(idx);
        end

    end
    methods %inherited
        function st=getSpikeTimes(obj) 
            st=obj.SpikeTable.SpikeTimes;
        end
        function sc=getSpikeClusters(obj)
            sc=obj.SpikeTable.SpikeCluster;
        end
    end
end

