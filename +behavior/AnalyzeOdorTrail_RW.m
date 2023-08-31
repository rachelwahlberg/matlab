function [dphi,nearX,nearY] = AnalyzeOdorTrail_RW(x0,y0) 
% zidphi
% Nothing actually to do with odor trails, just using to get zidphi scores
% Adapted from Papale's AnalyzeOdorTrail code
% 2017-07-12 AndyP
% Analyze one session 
% 2022-03-02 RachelW

C = [];
idphi = [];
dtime = 1/50;
window = 1;
postSmoothing = 0.5;

%[x0,y0,nx0,ny0,~,~] = Process_VT(position_results,startFrame);

%% (he didn't actually use) for each nan, find nearest position to each nan...
%// Index array for factor
x1 = 1:numel(x0);
%// Indices of NaNs
t2 = find(~isnan(x0));
%// Replace NaNs with the closest non-NaNs 
nearX0 = interp1(x1(t2),x0(t2),x1,'nearest');
nearY0 = interp1(x1(t2),y0(t2),x1,'nearest');

nearX0(1)=nearX0(2);
nearY0(1)=nearY0(2);

nearX = nearX0';% = cat(1,nearX,nearX0');
nearY = nearY0';% = cat(1,nearY,nearY0');
%% Get tortuosity and zidphi
dx = dxdt(nearX,dtime,window,postSmoothing); %x0 originally, confirm that it's working
dy = dxdt(nearY,dtime,window,postSmoothing);
%V0 = sqrt(dx.^2+dy.^2)./11.2;
%V = V0;% cat(1,V,V0); % cm/s

C0 = Tortuosity(dx,dy,dtime,window,postSmoothing);
C = cat(1,C,C0);
zC = nanzscore(abs(C));

% not yet integrated or z scored, despite function name. gives angular
% velocity
dphi = zIdPhi(dx,dy,dtime,window,postSmoothing); %dphi0

% Should be outside of the function
%idphi = cat(1,idphi,idphi0);
%zidphi=nanzscore(abs(idphi));




