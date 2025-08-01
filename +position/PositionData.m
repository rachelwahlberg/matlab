classdef PositionData < neuro.basic.ChannelTimeData & ...
        matlab.mixin.indexing.RedefinesParen
    %POSITIONDATA Summary of this class goes here
    %   Detailed explanation goes here

    properties
        units
        Info
    end

    methods
        function obj = PositionData(X,Y,Z,time)
            %PositionData Construct an instance of this class
            %   ticd should be in TimeIntervalCombined foormat
            if nargin>0
                % If the input X is already an instance of PositionData, copy its properties
                if isa(X, 'optiTrack.PositionData')
                    obj.data = X.data;
                    obj.time = X.timeIntervalCombined;
                    obj.units = X.units;
                    return;
                end

                % Validate the inputs
                if numel(X) ~= numel(Y) || numel(Z) ~= numel(Y) || numel(X) ~= time.getNumberOfPoints()
                    error('Sizes of XYZ or time are not equal.');
                end

                % Set the data property
                data = [X(:), Y(:), Z(:)];
                obj.data = array2table(data, 'VariableNames', {'X', 'Y', 'Z'});

                % Set the time property
                obj.time = time;

                % Set the units property
                obj.units = 'cm';
            end
        end

        function theta = getPositionDataAngle(obj)
            % for circular track - angle of data relative to center pt
            %data has come in 100x100, so move 0,0 to be center
            x = obj.data.X - 50;
            z = obj.data.Z - 50;

            theta = atan2(z,x);
            theta(theta < 0) = theta(theta < 0) + 2*pi; %have all between 0 and 2pi
            % figure %to make sure looks correct.
            % plot(theta)
        end

        function [positionData, idx] = getWindow(obj, range)
            % GETWINDOW Returns a new PositionData object and its corresponding index
            %   for a specified time range.
            %
            %   Inputs:
            %       obj         - a PositionData object
            %       range       - a 2-element vector specifying the start and end times
            %                     of the desired time range
            %
            %   Outputs:
            %       positionData - a new PositionData object containing the data within the
            %                      specified time range
            %       idx          - the index of the data points within the specified time range

            % Get the TimeIntervalData object associated with the PositionData object
            ticd = obj.time;

            % Create a copy of the input PositionData object
            positionData = obj;

            % Get the time window corresponding to the specified range
            window = ticd.getTimeIntervalForTimes(range);

            % Update the time property of the new PositionData object
            positionData.time = window;

            % Get the sample indices corresponding to the specified range
            samples = ticd.getSampleForClosest(range);

            % Get the data points within the specified range
            idx = samples(1):samples(2);
            positionData.data = obj.data(idx,:);
        end
        function obj=plus(obj,pd)
            if obj.time.getEndTime>=pd.time.getStartTime
                error(['Position data starts(%s) before the ' ...
                    'original ends(%s).'],pd.time.getStartTime, ...
                    obj.time.getEndTime)
            end
            if obj.time.getSampleRate~=pd.time.getSampleRate
                error('Samplerates are different. First:%dHz, Second:%dHz', ...
                    obj.time.getSampleRate,pd.time.getSampleRate)
            end
            newtime=obj.time+pd.time;
            newdata=[obj.data; pd.data];
            obj.time=newtime;
            obj.data=newdata;
        end
        function data=getData(obj)
            data=obj.data;
        end
        function obj=setData(obj,data)
            obj.data=data;
        end
        function mat=flatten2(obj)
            mat=obj.data(1:2,:);
        end
        function mat=flatten3(obj)
            mat=obj.data(1:3,:);
        end

        function pdman = getManifold(obj)
            time = obj.time;
            timestr = matlab.lang.makeValidName(time.tostring);
            file1 = java.io.File(obj.source);
            manifoldFile = fullfile(char(file1.getParent), ...
                ['position.PositionDataManifold' timestr '.mat']);

            if exist(manifoldFile, 'file')
                s = load(manifoldFile);
                pdman = s.pdman;
                if isempty(pdman.units)
                    pdman.units = 'cm';
                end
                return;
            else
                manifoldFile1 = fullfile(char(file1.getParent), ...
                ['position.PositionDataManifold' '*' '.mat']);
                fs=dir(manifoldFile1);
                if numel(fs)==1
                    s = load(fullfile(fs.folder,fs.name));
                    pdman = s.pdman;
                    if isempty(pdman.units)
                        pdman.units = 'cm';
                    end
                    return;
                end
            end

            try
                close(123);
            catch
            end
            figure(123);
            f = gcf;
            f.Position(3:4) = [2500 1500];
            tiledlayout(2, 3);

            c.numberOfPoints = 300;
            c.neighbors = 7;

            sampleRate = obj.time.getSampleRate;
            obj1 = obj.getDownsampled(sampleRate * 0.1);
            spd = obj1.getSpeed(3).Values;
            [~, I] = sort(spd, 'descend', 'MissingPlacement', 'last');

            data1 = table2array(obj1.data)';
            data2 = data1(:, I(1:(size(data1, 2) / 20)));
            data2(:, any(isnan(data2))) = [];

            nexttile(1, [2 1]);
            obj.plot3DtimeContinuous;
            title('Original');

            nexttile(3, [1 1]);
            manifold = external.Manifold.Manifold("Description");
            manifold = manifold.createGraph(data2, 'verbose', 'neighbors', c.neighbors, 'numPoints', c.numberOfPoints);
            manifold.plotGraph;
            ax = gca;
            ax.DataAspectRatio = [1 1 1];
            title('Graph');

            manifold = manifold.shortestPath('verbose');

            nexttile(6, [1 1]);
            manifold = manifold.scale('plot', 'sammon');
            ax = gca;
            ax.DataAspectRatio = [1 1 1];
            title('Scaled');

            pdman = position.PositionDataManifold(obj, manifold);
            pdman.config = c;

            nexttile(2, [2 1]);
            pdman.plot3DtimeContinuous;
            title('Dimension Reduced');

            ff = logistics.FigureFactory.instance(char(file1.getParent));
            ff.save(['position-PositionDataManifold-' timestr]);
            save(manifoldFile, 'pdman', '-mat');
        end
        function [velocity] = getSpeed(obj, smoothingWindowInSeconds)
            data = table2array(obj.getData)'; %get x,y,z coordinates
            timeDiffSeconds = diff(seconds(obj.time.getTimePoints));
            timeDiffSeconds2=[timeDiffSeconds median(timeDiffSeconds)];
            squaredDiffs = zeros(size(data));
            for dimIndex = 1:size(data, 1)
                squaredDiffs(dimIndex, 1:(end-1) )= diff(data(dimIndex, :)).^2;
            end %combines over x,y,z dimensions
            speeds = sqrt(sum(squaredDiffs, 1))./timeDiffSeconds2;
            if exist('smoothingWindowInSeconds', 'var')
                speeds = smoothdata(speeds, 'gaussian', obj.time.getSampleRate * ...
                    smoothingWindowInSeconds);
            end
            velocity = neuro.basic.Channel('Velocity', speeds, obj.time);
        end
        function [velocity] = getSpeed1D(obj, smoothingWindowInSeconds)
            data = table2array(obj.getData)'; %get x,y,z coordinates
            data(2,:) = []; %get rid of y 
            timeDiffSeconds = diff(seconds(obj.time.getTimePoints));
            timeDiffSeconds2=[timeDiffSeconds median(timeDiffSeconds)];
            squaredDiffs = zeros(size(data));
            for dimIndex = 1:size(data, 1)
                squaredDiffs(dimIndex, 1:(end-1) )= diff(data(dimIndex, :)).^2;
            end %combines over x,y,z dimensions
            speeds = sqrt(sum(squaredDiffs, 1))./timeDiffSeconds2;
            if exist('smoothingWindowInSeconds', 'var')
                speeds = smoothdata(speeds, 'gaussian', obj.time.getSampleRate * ...
                    smoothingWindowInSeconds);
            end
            velocity = neuro.basic.Channel('Velocity', speeds, obj.time);
        end
        function [om]= getOccupancyMap(obj,xedges,zedges)
            if nargin==1
                om=neuro.placeField.OccupancyMap(obj,obj.time.getSampleRate);
            else
                om=neuro.placeField.OccupancyMap(obj, ...
                    obj.time.getSampleRate,xedges,zedges);
            end
            om.Units=obj.units;
        end
        function ax = plot(obj)
            numPointsInPlot = 100000;
            time = obj.time;
            t_org = seconds(time.getTimePoints() - (time.getZeitgeberTime() - time.getStartTime()));
            downsampleFactor = max(1, round(numel(t_org) / numPointsInPlot)); % add lower bound
            data = table2array(obj.getData());
            downsampledData = downsample(data, downsampleFactor); % remove loop
            t = hours(seconds(downsample(t_org, downsampleFactor)));
            plot(t, downsampledData);
            legend(obj.getData().Properties.VariableNames);
            xlabel('ZT (Hrs)');
            ylabel(['Location (', obj.units, ')']);
            ax = gca;
        end

        function ax = plot2D(obj, numPointsInPlot)
            % Set default value for numPointsInPlot if not provided
            if nargin < 2
                numPointsInPlot = 10000;
            end

            % Get time data
            ticd = obj.time;
            t_org = ticd.getTimePointsZT;

            % Calculate downsample factor
            downsampleFactor = round(numel(t_org) / numPointsInPlot);

            % Get input data
            dims = {'X','Z'};
            inputData = table2array(obj.data(:, dims));

            % Downsample input data
            downsampledData = downsample(inputData, downsampleFactor)';

            % Get colormap
            colorMap = linspecer(size(downsampledData, 2));

            % Plot scatter plot
            scatter(downsampledData(1,:), downsampledData(2,:), [], colorMap, ...
                'filled', 'MarkerFaceAlpha', .2, 'MarkerEdgeAlpha', .2, 'SizeData', 5);

            % Set axis labels and aspect ratio
            ylabel([dims{2} ' ' obj.units]);
            xlabel([dims{1} ' ' obj.units]);
            ax = gca;
            ax.DataAspectRatio = [1 1 1];

            % Set colormap
            colormap(colorMap);
        end
        function p = plot2DContinuous(obj, numPointsInPlot)
            if ~exist('numPointsInPlot','var')
                numPointsInPlot=10000;
            end

            % Get data and time information
            data = obj.getData();
            time = obj.time;
            t_org = time.getTimePointsZT();

            % Downsample data
            downsampleFactor = round(numel(t_org) / numPointsInPlot);
            if downsampleFactor < 1
                downsampleFactor = 1;
            end
            data = downsample(data, downsampleFactor);
            t = downsample(t_org, downsampleFactor);

            % Plot continuous lines
            dims = {'X','Z'};
            p = plot(data.(dims{1}), data.(dims{2}));
            p.LineWidth = 2;
            p.Color = 'k';
            ylabel([dims{2} ' ' obj.units]);
            xlabel([dims{1} ' ' obj.units]);
        end
        function ax = plot3Dtime(obj, numPointsInPlot)
            % Set default value for numPointsInPlot if not provided
            if nargin < 2
                numPointsInPlot = 10000;
            end

            % Get time vector and downsample factor
            ticd = obj.time;
            t_org = seconds(ticd.getTimePointsZT);
            downsampleFactor = max(round(numel(t_org) / numPointsInPlot), 1);

            % Get data and downsample along the first two dimensions
            dims = {'X', 'Z', 'Y'};
            data0 = obj.getData();
            data1 = table2array(data0(:, dims))';
            data2 = zeros(size(data1));
            for ich = 1:size(data1, 1)
                data2(ich, :) = downsample(medfilt1(data1(ich,:), ...
                    ticd.getSampleRate()), downsampleFactor);
            end

            % Add time to the third dimension
            data2(3, :) = data2(3, :) + linspace(1, t_org(end) - t_org(2), size(data2, 2));

            % Plot data as a 3D scatter plot
            color1 = linspecer(size(data2, 2));
            scatter3(data2(1,:), data2(2,:), data2(3,:), [], color1, ...
                'filled', 'MarkerFaceAlpha', .2, 'MarkerEdgeAlpha', .2, 'SizeData', 5);
            xlabel(dims{1});
            ylabel(dims{2});
            zlabel('Time (s)');
            ax = gca;
            ax.DataAspectRatio = [1 1 4];
        end

        function ax = plot3DtimeContinuous(obj, numPointsInPlot)
            if ~exist('numPointsInPlot','var')
                numPointsInPlot=10000;
            end
            colors=colororder;
            ticd=obj.time;
            t_org=seconds(ticd.getTimePointsZT);
            downsamplefactor=round(numel(t_org)/numPointsInPlot);
            if downsamplefactor<1
                downsamplefactor=1;
            end

            dims={'X','Z','Y'};
            data0=obj.getData;
            data1=table2array(data0(:,dims))';
            for ich=1:size(data1,1)
                data2(ich,:)=downsample(data1(ich,:),downsamplefactor); %#ok<AGROW> %every 5th sample
            end
            lintime=linspace(0,t_org(end)-t_org(2),size(data2,2));
            data2(3,:)=data2(3,:)+lintime; % add time
            plot3(data2(1,:),data2(2,:),data2(3,:),LineWidth=1,Color= ...
                colors(1,:));
            hold on
            zlabel('Time (s)');
            ylabel(dims{2});
            xlabel(dims{1});
            ax=gca;
            ax.DataAspectRatio=[1 1 4];
            ax.ZDir="reverse";
            ax.YDir="reverse";
            %             nanidx=any(isnan(data2));
            %             nanarr=nan(size(data2));
            %             nanarr(1:2,nanidx)=0;
            %             nanarr(3,nanidx)=lintime(nanidx);
            %             plot3(nanarr(1,:),nanarr(2,:),nanarr(3,:),
            % LineWidth=2,Color=colors(2,:));
            %             data2(1:2,:)=0;
            %             data2(3,~nanidx)=lintime(~nanidx);
            %             plot3(data2(1,:),data2(2,:),data2(3,:),
            % LineWidth=2,Color=colors(3,:));
        end
        function ax = plot3D(obj, numPointsInPlot)
            if ~exist('numPointsInPlot','var')
                numPointsInPlot=10000;
            end
            ticd=obj.time;
            t_org=ticd.getTimePointsZT;
            downsamplefactor=round(numel(t_org)/numPointsInPlot);

            dims={'X','Z','Y'};
            data0=obj.getData;
            data1=table2array(data0(:,dims))';
            for ich=1:size(data1,1)
                data2(ich,:)=downsample(medfilt1(data1(ich,:), ...
                    ticd.getSampleRate),downsamplefactor); %#ok<AGROW>
            end
            plot3(data2(1,:),data2(2,:),data2(3,:));
            zlabel(dims{3});
            ylabel(dims{2});
            xlabel(dims{1});
        end
        function [] = plot3DMark(obj,event,color)
            if ~exist('color','var')
                color=[];
            end
            ticd=obj.time;
            t_org=seconds(ticd.getTimePoints);

            dims={'X','Z','Y'};
            data0=obj.getData;
            data1=table2array(data0(:,dims))';
            data1(3,:)=data1(3,:)+linspace(1,t_org(end)-t_org(2), ...
                size(data1,2)); % add time

            data2=data1(:,event);
            try
                s=scatter3(data2(1,:),data2(2,:),data2(3,:),[],color,'filled', ...
                    'MarkerFaceAlpha',.5,'MarkerEdgeAlpha',.5,'SizeData',30);
            catch ME
                s=scatter3(data2(1,:),data2(2,:),data2(3,:),'filled', ...
                    'MarkerFaceAlpha',.5,'MarkerEdgeAlpha',.5,'SizeData',30);
            end
                 if isempty(color)
                s.MarkerFaceColor="k";
            end
            ylabel([dims{2} ' ' obj.units]);
            xlabel([dims{1} ' ' obj.units]);

        end

        function [] = plot3DEventPlane(obj,event,color)
            % event comes in as indices (see SpikeUnitArray.plot3Dmark for
            % ex)
             if ~exist('color','var')
                color=[0.3010 0.7450 0.9330];
            end
            ticd=obj.time;
            t_org=seconds(ticd.getTimePoints);

            dims={'X','Z','Y'};
            data0=obj.getData;
            data1=table2array(data0(:,dims))';
            data1(3,:)=data1(3,:)+linspace(1,t_org(end)-t_org(2), ...
                size(data1,2)); % add time

            data2=data1(:,event);

            xvec = linspace(min(data1(1,:)),max(data1(1,:))); % in X space
            zvec = linspace(min(data1(2,:)),max(data1(2,:))); % in Y space
            
           plane = ones(length(xvec),length(zvec))*data2(3);
           m = mesh(xvec,zvec,plane,'FaceAlpha',0.2,'EdgeAlpha',0.3,'EdgeColor',color,'FaceColor',color);
        end

        function [] = plot2DMark(obj,eventtomark) %mark refers to an event such as a timestamp
            dims={'X','Z','Y'};
            data0=obj.getData;
            data1=table2array(data0(:,dims))';

            data2=data1(:,eventtomark);
            s=scatter(data2(1,:),data2(2,:),[],'filled', ...
                'MarkerFaceAlpha',.5,'MarkerEdgeAlpha',.5,'SizeData',5);
            s.MarkerFaceColor="k";
            ylabel([dims{2} ' ' obj.units]);
            xlabel([dims{1} ' ' obj.units]);

        end
        function obj = getTimeWindow(obj,timeWindow)
            %can be relative, zt, or absolute time object (see
            %time.TimeIntervalAbstract)
            ticd=obj.time;
            ticdnew=ticd.getTimeIntervalForTimes(timeWindow);
            s1=ticd.getSampleForClosest(ticdnew.getStartTime);
            s2=ticd.getSampleForClosest(ticdnew.getEndTime);
            obj.time=ticdnew;
            obj.data=obj.data(s1:s2,:);
        end
        function obj = getDownsampled(obj,dsfactor)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            data1=table2array(obj.data)';
            for idim=1:size(data1,1)
                data1ds(idim,:)=downsample(medfilt1(data1(idim,:),dsfactor), ...
                    dsfactor,dsfactor-1);
            end
            obj.data=array2table(data1ds',"VariableNames", ...
                obj.data.Properties.VariableNames);
            ticd=obj.time;
            obj.time=ticd.getDownsampled(dsfactor);
        end
        function obj = getMedianFiltered(obj,winseconds)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            data1=table2array(obj.data)';
            data2=smoothdata(data1,2,'movmedian',winseconds*obj.time.getSampleRate);
            obj.data=array2table(data2',"VariableNames", ...
                obj.data.Properties.VariableNames);
        end
        function obj = getMeanFiltered(obj,winseconds)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            data1=table2array(obj.data)';
            data2=smoothdata(data1,2,'movmean',winseconds*obj.time.getSampleRate);
            obj.data=array2table(data2',"VariableNames", ...
                obj.data.Properties.VariableNames);
        end
        function [data1,sample] = getPositionForTimes(obj,times)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            if numel(times)>0
                time=obj.time;
                sample=time.getSampleForClosest(times);
                sample = sample(~isnan(sample));
                data1=obj.data(sample,:);
            else
                sample=[];
                data1=obj.data(sample,:);
            end
        end

        function [positions] = getPositions_largeArrays(obj,ztTimes,stSampleRate)
            %A more efficient way to get positions for large arrays of
            %spiketimes compared to getPositionForTimes.

            %pos=obj.PositionData;
            pt=seconds(obj.time.getTimePointsZT)'; % get positionData timepoints ZT
            sts=seconds(ztTimes); %all spikes, not divided by label
            % Create a KDTree object from pt
            tree = KDTreeSearcher(pt);
            % Find the indices of the closest values in pt for each value
            % in sts % could be way off because a lot of it isn't on the
            % track itself? 
            locinpos = knnsearch(tree, sts);
            threshold=max(1./[obj.time.getSampleRate ...
                stSampleRate]); %max between 1/behavior sr and 1/lfp sr
            locinpos(abs(sts - pt(locinpos)) > threshold) = NaN;
            positions=array2table(nan(length(sts),3), ...
                VariableNames=obj.data.Properties.VariableNames);
            idxvalid=~isnan(locinpos);
            positions(idxvalid,:)=obj.data(locinpos(idxvalid),:);
        end

        function data1 = getPositionForSamples(obj,samples)
            data1=obj.data(samples,:);
        end

        function [obj, folder]= saveInPlainFormat(obj,folder,ext1)
            if ~exist('ext1','var')
                ext1='position.points.csv';
            end
            extt='position.time.csv';
            if exist('folder','var')
                if ~isfolder(folder)
                    folder= pwd;
                    warning(['The folder does not exist. ' ...
                        'Files will be created in %s.'],folder);
                end
            else
                folder= fileparts(obj.source);
            end
            time=obj.time; %#ok<*PROPLC>
            timestr=matlab.lang.makeValidName(time.tostring);
            time.saveTable(fullfile(folder,[timestr extt]));
            file1=fullfile(folder,[timestr ext1]);
            writetable(obj.data,file1);
            file2=fullfile(folder,[timestr '.mat']);
            save(file2,'obj');
            folder=string(py.os.path.realpath(py.os.path.expanduser(folder)));
            %             obj=obj.loadPlainFormat(folder);
        end
        function obj= loadPlainFormat(obj,folder)
            ext1='position.points.csv';
            extt='position.time.csv';
            [file1, uni]=obj.getFile(folder,ext1);
            obj.source=file1;
            obj.data=readtable(obj.source);
            folder=fileparts(file1);
            obj.time=time.TimeIntervalCombined( ...
                fullfile(folder,[uni extt]));
        end
        function [file2, uni]=getFile(~,folder,extension)
            if ~exist('folder','var')
                folder= pwd;
            end
            if isfile(folder)
                [folder1,name,ext1]=fileparts(folder);
                uni1=split([name ext1],extension);
                uni=uni1{1};
                file1=dir(fullfile(folder1,[uni,extension]));
            else
                file1=dir(fullfile(folder,['*' extension]));
                if numel(file1)>1
                    [name,folder1] = uigetfile({['*' extension],extension}, ...
                        'Selectone of the position files',folder);
                    file1=dir(fullfile(folder1,name));
                end
            end
            file2=fullfile(file1.folder,file1.name);
            [~,name,ext1]=fileparts(file2);
            uni1=split([name ext1],extension);
            uni=uni1{1};
        end

        function obj = normalizePositionData(obj)

        x = obj.data.X; z = obj.data.Z; y = obj.data.Y; %vertical dimension

        min_x = min(x); min_z = min(z); min_y = min(y);
        x = x - min_x; z = z - min_z; y = y-min_y; %shift to start at 0;

        max_x = max(x); max_z = max(z); max_y = max(y);
        x = x/max_x; z = z/max_z; y = y/max_y; %normalize to 1

        obj.data.X = x; obj.data.Z = z; obj.data.Y = y; %resave

        end


        function linearized = linearizePositionDataTemplate(obj,template)
            %template is currently manually made - put in template.X and
            %template.Y vectors for the template you would like to
            %linearize to. Spits out positionData object but with lin. data
            npoints = length(obj.data.X);

            x = obj.data.X;
            z = obj.data.Z;

            linearized = obj;
            %linearized.data.Y = zeros(size(linearized.data.Y)); %Y is irrelevant here
            for p = 1:npoints
                if isnan(x(p))
                    continue
                end

                d = sqrt((template.X(:)-x(p)).^2+(template.Z(:)-z(p)).^2);
                closest_idx = find(abs(d)==min(abs(d)));

                if isempty(closest_idx) %for noninterpolated data
                linearized.data.X(p) = template.X(last_closest);
                linearized.data.Z(p) = template.Z(last_closest);
                continue
                else
                linearized.data.X(p)=template.X(closest_idx(1));
                linearized.data.Z(p)=template.Z(closest_idx(1));
                end
                last_closest = closest_idx(1);
            end

            % figure
            % hold on
            % plot(positionData.data.X,positionData.data.Z,'.');
            % h=circle(circleParams.circCenter,circleParams.trackRadius,2000,'r'); %make a circle template
            % template.X=h.XData;
            % template.Z=h.YData;
        end

        function obj = eliminateOffCircTrack(obj,params,plotfig) %vals are specific to 6/14
            % specific for use with RW circular track - remaining points
            % will only be those on the track itself, getting rid of
            % sticking head over side/etc.

            %innerRadius = 0.35;
            %outerRadius = 0.41;
            %circCenter = [0.48,0.51];

            x = obj.data.X; z = obj.data.Z;            
           % npoints = 2000;

           if exist("plotfig","var")
            figure
            hold on
            plot(x,z,'.b');
            circle(params.circCenter,params.outerRadius,2000,'r'); %outer circle
            circle(params.circCenter,params.innerRadius,2000,'g'); %inner circle
           end

            %eliminate those within the inner circle and outside outer circle
            
            distFromCenter = sqrt((params.circCenter(1) - x).^2 + (params.circCenter(2) - z).^2);
            outRemove = distFromCenter > params.outerRadius;
            inRemove = distFromCenter < params.innerRadius;

            idxall = outRemove|inRemove;
            obj.data.X(idxall)=nan;
            obj.data.Y(idxall)=nan;
            obj.data.Z(idxall)=nan;
        end
    end % methods end
    methods (Access=protected)
        function obj = parenReference(obj, indexOp)
            try
                idx=indexOp.Indices{:};
            catch ME
                idx=indexOp(1,1).Indices{:};
            end
            data=table2array(obj.data);
            data(~idx,:) =nan;
            obj.data=array2table(data,VariableNames= ...
                obj.data.Properties.VariableNames);
        end

        function obj = parenAssign(obj,indexOp,varargin)
            % Ensure object instance is the first argument of call.
            if isempty(obj)
                obj = varargin{1};
            end
            if isscalar(indexOp)
                assert(nargin==3);
                rhs = varargin{1};
                obj.ContainedArray.(indexOp) = rhs.ContainedArray;
                return;
            end
            [obj.(indexOp(2:end))] = varargin{:};
        end

        function n = parenListLength(obj,indexOp,ctx)
            if numel(indexOp) <= 2
                n = 1;
                return;
            end
            containedObj = obj.(indexOp(1:2));
            n = listLength(containedObj,indexOp(3:end),ctx);
        end

        function obj = parenDelete(obj,indexOp)
            obj.ContainedArray.(indexOp) = [];
        end
    end
    methods
        function out = cat(dim,varargin)
            numCatArrays = nargin-1;
            newArgs = cell(numCatArrays,1);
            for ix = 1:numCatArrays
                if isa(varargin{ix},'ArrayWithLabel')
                    newArgs{ix} = varargin{ix}.ContainedArray;
                else
                    newArgs{ix} = varargin{ix};
                end
            end
            out = ArrayWithLabel(cat(dim,newArgs{:}));
        end

        function varargout = size(obj,varargin)
            [varargout{1:nargout}] = [1 1];
        end
    end
end

