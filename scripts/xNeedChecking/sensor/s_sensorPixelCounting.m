%% s_pixelCounting
%
% There is a relationship between pixel density and sensor size.  We
% calculate examples here using some standard sensor format sizes.
% Manufacturers can alter this a bit by using slightly different imaging
% areas.
%
% Copyright ImagEval Consultants, LLC, 2010

% You can cheat a little on the imaging area, but not the overall form
% factor 
cheatArea = 1;  

sensorSizeMeters = sensorFormats('quarterinch')*cheatArea;
pixelSizeMeters  = (1:.2:2.4)*1e-6;
N = length(pixelSizeMeters);

pixelCount = zeros(size(pixelSizeMeters));
for ii=1:N
    pixelCount(ii) = prod(sensorSizeMeters)/pixelSizeMeters(ii)^2;
    fprintf('Size: %f, Megapixels = %f\n',pixelSizeMeters(ii)*10^6, pixelCount(ii)/10^6);
end

%% End
