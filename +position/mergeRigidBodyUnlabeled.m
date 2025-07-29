

clear 
%%
%folder=['/data1/EphysAnalysis/SleepDeprivationData/' ...
   % 'AE_2019-10-23_NSD/_Position'];

folder='/data/ExperimentsRSW/CircularMaze/Harry/CA1/20230614/Optitrack';
olrb=position.optiTrack.OptiLoader.instance(folder);
olother=position.optiTrack.OptiLoader.instance(folder);
ofrb=olrb.loadFile;
ofother=olother.loadFile;
positionsrb=ofrb.Positions;
positionsother=ofother.Positions;
rbln=positionsrb(ismember(positionsrb.Type,'Rigid Body'),:);
rb=rbln.PositionData{1}; % get rigid body 

%% detect rigidbody nans
rbmlns=positionsother(ismember(positionsother.Type,'Marker') & ...
    contains(positionsother.Name,'Marker')...
    ,:);
clear nan1;
for irbm=1:height(rbmlns)
    rbmpd=rbmlns(irbm,"PositionData").PositionData{1};
    if exist('nan1','var')
        nan1=nan1 & isnan(rbmpd.data.X);
    else
        nan1=isnan(rbmpd.data.X);
    end
end
%% unlabeled marker estimations
umlns=positionsother(ismember(positionsother.Type,'Marker') & ...
    contains(positionsother.Name,'Unlabeled')...
    ,:);
clear matx maty matz;
for irbm=1:height(umlns)
    umpd=umlns(irbm,"PositionData").PositionData{1};
    if exist('matx','var')
        matx=[matx umpd.data.X];
        maty=[maty umpd.data.Y];
        matz=[matz umpd.data.Z];
    else
        matx=umpd.data.X;
        maty=umpd.data.Y;
        matz=umpd.data.Z;
    end
end
data1(:,1)=median(matx,2,"omitnan");
data1(:,2)=median(maty,2,"omitnan");
data1(:,3)=median(matz,2,"omitnan");
datatbl=array2table(data1,VariableNames={'X','Y','Z'});
%% rigid body missing points replacement

pos1=rb;
pos1.data(nan1,:)=array2table(nan(length(find(nan1)),3));
pos1.data(nan1,:)=datatbl(nan1,:);
pos2=pos1;

% adjust to the actual size of the track
%scatter(pos2.data.X,pos2.data.Z)% to see top view
%scatter(pos2.data.X,pos2.data.Y)%to see side view

pos2.data.X=pos2.data.X    *23;
pos2.data.Y=pos2.data.Y    *23;
pos2.data.Z=pos2.data.Z    *23;
%idx1=pos2.data.X<-54|pos2.data.X>56;
%idx2=pos2.data.Y<45|pos2.data.Y>65;
%idx3=pos2.data.Z<-93|pos2.data.Z>90;
idx1=pos2.data.X<-60|pos2.data.X>73;
idx2=pos2.data.Y<-20|pos2.data.Y>34;
idx3=pos2.data.Z<-10|pos2.data.Z>121;
idxall=idx1|idx2|idx3;
pos2.data.X(idxall)=nan;
pos2.data.Y(idxall)=nan;
pos2.data.Z(idxall)=nan;
pos3=pos2;
pos3.data=fillmissing(pos3.data,"makima",MaxGap=30*pos3.time.getSampleRate,EndValues="nearest");
pos3=pos3.getMedianFiltered(.5);
pos=pos3.getLowpassFiltered(6);
pos.saveInPlainFormat(folder)