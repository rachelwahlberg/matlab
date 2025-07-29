function [] = distanceontrack(angle,tracklength)

%angle in degrees

distanceInInches = angle/180*(tracklength*12);
distanceRounded = round(distanceInInches/(1/32))/32;

inchesRounded = floor(distanceRounded);
inchesFraction = (distanceRounded - inchesRounded)/(1/32);

disp([num2str(inchesRounded) ' inches and ' num2str(inchesFraction) '/32'])

