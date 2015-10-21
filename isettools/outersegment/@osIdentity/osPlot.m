function osPlot(obj, sensor)
% osPlot: a method of @osIdentity to plot the computed results of the 
% outsergment object.
%  
% Inputs: the osIdentity object and the sensor object.
% 
% Outputs: no variables, but a figure with two subplots is generated. The
% first shows the input signal and the second shows the output signal.
% 
% osPlot(identityOS, sensor);
% 
% 8/2015 JRG NC DHB

% fprintf('<strong>\n%s:\n\t%s()\n</strong>', class(obj), mfilename());

dt = sensorGet(sensor, 'time interval');

% Plot input signal (isomerizations) at a particular (x, y) over time.
figNum = 1; 
h = vcNewGraphWin([],'wide');
set(h, 'Name', sprintf('Output of %s', class(obj)));

% Plot output signal at a particular (x, y) over time.

[sz1 sz2 sz3] = size(obj.rgbData); 
outputSignal = squeeze(obj.rgbData(round(sz1/2),round(sz2/2),:,:));
plot((0:numel(outputSignal(:,1))-1)*dt, outputSignal, 'k-');

title('output signal, RGB');
xlabel('Time (sec)');
ylabel('Luminance');
drawnow;