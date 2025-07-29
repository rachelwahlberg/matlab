function positionData = fillRigidBodyNans(pos0,savefolder)
% 
% folder='/data/ExperimentsRSW/CircularMaze/Harry/CA1/20230614/Optitrack';
% 
% olrb=position.optiTrack.OptiLoader.instance(folder);
% ofrb=olrb.loadFile;
% positionsrb=ofrb.Positions;
% 
% olother=position.optiTrack.OptiLoader.instance(folder);
% ofother=olother.loadFile;
% pos0=ofother.Positions;

%% detect rigidbody nans
rbmlns=pos0(ismember(pos0.Type,'Marker') & ...
 contains(pos0.Name,'Marker'),:); %get other markers
% rbmlns= pos0(contains(pos0.Type,'Rigid') & ...
%     contains(pos0.Name,'Marker'),:); %for Rigid Body

for irbm=1:height(rbmlns) %get the nans from all the rigid bodies
    rbmpd=rbmlns(irbm,"PositionData").PositionData{1};
    if exist('nan1','var')
        nan1=nan1 & isnan(rbmpd.data.X);
    else
        nan1=isnan(rbmpd.data.X);
    end
end
%% unlabeled marker estimations % combines all the unlabled markers, takes the median of them.
umlns=pos0(contains(pos0.Type,'Marker') & ...
    contains(pos0.Name,'Unlabeled'),:);

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
rbln=pos0(ismember(pos0.Type,'Rigid Body'),:);
%rbln=pos0(contains(pos0.Name,'RigidBody'),:); % UK DATA
if height(rbln)>1
    for i = 1:height(rbln)
        allX(:,i) = rbln.PositionData{i}.data.X;
        allY(:,i) = rbln.PositionData{i}.data.Y;
        allZ(:,i) = rbln.PositionData{i}.data.Z;
    end

    rb=rbln.PositionData{1}; % get rigid body
    rb.data.X = mean(allX,2,"omitnan");
    rb.data.Y = mean(allY,2,"omitnan");
    rb.data.Z = mean(allZ,2,"omitnan");
else
    rb=rbln.PositionData{1};
end

pos1=rb;
pos1.data(nan1,:)=array2table(nan(length(find(nan1)),3)); %make all the nans consistent across rows I think
pos1.data(nan1,:)=datatbl(nan1,:);
pos2=pos1;

% adjust to the actual size of the track
%scatter(pos2.data.X,pos2.data.Z)% to see top view
%scatter(pos2.data.X,pos2.data.Y)%to see side view

%idxall =pos2.data.Z <0 & pos2.data.X > -200; % UK DATA
%%%%%%%%% bring back for RW data %%%%%
pos2.data.X=pos2.data.X    *23;
pos2.data.Y=pos2.data.Y    *23;
pos2.data.Z=pos2.data.Z    *23;
idx1=pos2.data.X<-60|pos2.data.X>73;
idx2=pos2.data.Y<-20|pos2.data.Y>34;
idx3=pos2.data.Z<-10|pos2.data.Z>121;
idxall=idx1|idx2|idx3;
%%%%%%%%%%%%%%%%%%%%%%%%%
pos2.data.X(idxall)=nan;
pos2.data.Y(idxall)=nan;
pos2.data.Z(idxall)=nan;
pos3=pos2;
pos3.data=fillmissing(pos3.data,"makima",MaxGap=30*pos3.time.getSampleRate,EndValues="nearest"); %if maxgap is too low, can result in nans getting through. 
pos3=pos3.getMedianFiltered(.5);
positionData=pos3.getLowpassFiltered(6);

if ~isempty(savefolder)
positionData.saveInPlainFormat(savefolder)
end


end