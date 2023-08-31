% VTE Analyses pipeline

%% Load data


%% Get indices for general vte region
vte_coord.x = [-14000 -8000];
vte_coord.y = [-2000 2000];

nsmps = length(optitrack.position.x);
x_ind = zeros(1,nsmps);
y_ind= zeros(1,nsmps);
x_ind = (optitrack.position.x > vte_coord.x(1)) .* (optitrack.position.x < vte_coord.x(2));
y_ind = (optitrack.position.y > vte_coord.y(1)) .* (optitrack.position.y < vte_coord.y(2));

region_ind = logical(x_ind .* y_ind);

x_pos = nan(size(optitrack.position.x));
y_pos = nan(size(optitrack.position.y));
x_pos(region_ind) = optitrack.position.x(region_ind);
y_pos(region_ind) = optitrack.position.y(region_ind);

x_ori = nan(size(optitrack.position.x));
y_ori = nan(size(optitrack.position.y));
x_ori(region_ind) = optitrack.orientation.x(region_ind);
y_ori(region_ind) = optitrack.orientation.y(region_ind);



figure
plot(x_pos,y_pos);

plot(optitrack.position.x(region_ind),optitrack.position.y(region_ind))

% next would be identifying what a vte position is in itself (rotation
% based?)
% A vte would be rotation in the y direction
% positive rotation than negative rotation and vice versa?


