function [uData,g] = plotSensorEtendue(sensor)
%Mesh representing etendue entries in the ISA structure
%
%   [uData,g] = plotSensorEtendue([sensor])
%
% The etendue is computed using mlAnalyzeArrayEtendue.
%
% Examples
%  plotSensorEtendue(vcGetObject('ISA'));
%  plotSensorEtendue;
%
% Copyright ImagEval Consultants, LLC, 2005.

if notDefined('sensor'), sensor = vcGetObject('sensor'); end

% Make a figure showing the etendue across the array.  The units of the
% support are unclear to me at this moment.
g = vcNewGraphWin;
set(g,'name',sprintf('ISET: Etendue (%s)',sensorGet(sensor,'vignettingName')));

sensorEtendue = sensorGet(sensor,'etendue');
if isempty(sensorEtendue), 
    sensor = sensorVignetting(sensor); 
    sensorEtendue = sensorGet(sensor,'etendue');
end

support = sensorGet(sensor,'spatialSupport','um');
mesh(support.x,support.y,sensorEtendue);
zRange = get(gca,'zlim');
if zRange(2) - zRange(1) < 0.1;
    zRange(1) = floor(zRange(1)*10)/10;
    zRange(2) = zRange(1) + 0.2;
    set(gca,'zlim',[zRange(1) zRange(2)]);
end
xlabel('Position (um)')
ylabel('Position (um)')
zlabel('Relative illumination')

uData.support = support;
uData.sensorEtendue = sensorEtendue;
set(g,'userdata',uData);

return;
