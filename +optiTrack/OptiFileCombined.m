classdef OptiFileCombined < time.Timelined
    %OPTIFILECOMBINED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        OptiFiles
    end
    
    methods
        function obj = OptiFileCombined(varargin)

            optiFiles=utilities.CellArrayList();


            if iscell(varargin{1}) % RW add
                files = varargin{1};
            else
                files = varargin;
            end

            for iArgIn=1:length(files)
                optifile=files{iArgIn};
                assert(isa(optifile,'optiTrack.OptiFile'));
                optiFiles.add(optifile);
            end

            obj.OptiFiles=optiFiles;
        end
        function obj=plus(obj,varargin)
            for iArgIn=1:(nargin-1)
                optifile=varargin{iArgIn};
                assert(isa(optifile,'optiTrack.OptiFile'));
                obj.OptiFiles.add(optifile);
            end
            obj=obj.getSorted;
        end
        function tls=getTimeline(obj)
            iter=obj.getIterator();
            tls=[];
            i=1;
            while(iter.hasNext)
                optifile=iter.next();
                tl=optifile.getTimeline();
                tls{i}=tl;i=i+1;
            end
        end
        function ofs=getOptiFiles(obj)
            ofs=obj.OptiFiles;
        end
        function positionData=getMergedPositionData(obj)
            ofs=obj.getSorted;
            for iof=1:ofs.OptiFiles.length
                of=ofs.OptiFiles.get(iof);
                if isa(of,'optiTrack.OptiCSVFileSingleMarker')
                    pd=of.getMergedPositionData;
                    of.table=pd.data;
                end
                ti1=of.getTimeInterval;
                if exist('table1','var')
                    table1=[table1;of.table(:,{'X','Y','Z'})];
                    tic1=tic1+ti1;
                else
                    table1=of.table(:,{'X','Y','Z'});
                    tic1=ti1;
                end
            end
            X=table1.X;
            Y=table1.Y;
            Z=table1.Z;
            prompt = {'Zeitgeber Time:'};
            dlgtitle = datestr(tic1.getDate);
            dims = [1 10];
            definput = {'08:00'};
            zt = duration(inputdlg(prompt,dlgtitle,dims,definput),'InputFormat','hh:mm');
            tic1=tic1.setZeitgeberTime(zt);
            positionData=optiTrack.PositionData(X,Y,Z,tic1);
        end

        function positionData = getInterpolatedMergedPosition(obj)
            positionData = getMergedPositionData(obj);
            x = positionData.data.X;
            y = positionData.data.Y;
            z = positionData.data.Z;

            x1 = 1:numel(x);
            x2 = find(~isnan(x)); %// Indices of NaNs

            %// Replace NaNs with the closest non-NaNs
            positionData.data.X = interp1(x1(x2),x(x2),x1,'nearest'); % interpolated values!
            positionData.data.Y = interp1(x1(x2),y(x2),x1,'nearest'); % interpolated values
            positionData.data.Z = interp1(x1(x2),z(x2),x1,'nearest');
    
        end

        function obj=getSorted(obj)
            ofs=obj.OptiFiles;
            for iof=1:ofs.length
                of=ofs.get(iof);
                list(iof)=datenum(of.getStartTime);
            end
            [B,I] = sort(list);
            ofssorted=utilities.CellArrayList();
            for iof=1:ofs.length
                ofssorted.add(ofs.get(I(iof)));
            end
            obj.OptiFiles=ofssorted;
        end
    end
    methods (Access=private)
        function iterator=getIterator(obj)
            iterator=obj.OptiFiles.createIterator;
        end
        
    end
end

