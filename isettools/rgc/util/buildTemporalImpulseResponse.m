function k = buildTemporalImpulseResponse(sensor, namesCellTypes)
% buildTemporalImpulseResponse: a util function of the @rgc parent class,
% this builds the temporal filter function.

% From code by J. Pillow:
% Create a default (temporal) stimulus filter

integrationTime = sensor.integrationTime;

% DTsim = .01; % Bin size for simulating model & computing likelihood.
% nkt = 20;  % Number of time bins in filter;
DTsim = integrationTime; % Bin size for simulating model & computing likelihood.
filterLength = 0.2;
nkt = ceil(filterLength/integrationTime);  % Number of time bins in filter;
timeAxis = integrationTime:integrationTime:filterLength;
tk = [0:nkt-1]';
b1 = nkt/32; b2 = nkt/16;
k1 = 1/(gamma(6)*b1)*(tk/b1).^5 .* exp(-tk/b1);  % Gamma pdfn
k2 = 1/(gamma(6)*b2)*(tk/b2).^5 .* exp(-tk/b2);  % Gamma pdf
k = flipud(k1-k2./1.5);
% figure; plot(timeAxis,k);