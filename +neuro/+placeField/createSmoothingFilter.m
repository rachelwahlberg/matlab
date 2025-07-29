function ff = createSmoothingFilter(rr) %r 1 as default

[xx,zz] = meshgrid(-2*rr:2*rr,-2*rr:2*rr);
ff = exp(-((xx).^2 + (zz).^2)/4); %the filter
ff = ff/(sum(ff(:)));

end