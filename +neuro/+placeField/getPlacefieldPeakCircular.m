function peak = getPlacefieldPeakCircular(placemap)
% R Wahlberg, based on Kaya code (in PlaceFieldMap class)
% for use with linearized circular track, 
% place map normalized between 0 and 1 but corresponding to 0 to 2pi rang

% outputs peak, a table with FiringRate and Position variables.

[pks,locs1]=max(placemap);
angle=linspace(0,2*pi,length(placemap));
locs=angle(locs1);
peak=table(pks',locs', ...
    VariableNames={'FiringRate','Position'});

end