function osBioPhysPlot(obj, sensor)
% osBioPhysPlot: a method of @osBioPhys to plot the computed results of the 
% outsergment object.
% 
% 
% 8/2015 JRG NC DHB

% fprintf('<strong>\n%s:\n\t%s()\n</strong>', class(obj), mfilename());

dt = sensorGet(sensor, 'time interval');

% Plot input signal (isomerizations) at a particular (x, y) over time.
figNum = 1; 
h = figure();
set(h, 'Name', sprintf('Output of %s', class(obj)));
set(h, 'Position', [10+50*figNum 10+50*figNum, 1024 256]);
subplot(1,2,1);
[sz1 sz2 sz3] = size(sensor.data.volts); 
inputSignal(1,:) = sensor.data.volts(round(sz1/2),round(sz2/2),:);
plot((0:numel(inputSignal)-1)*dt, inputSignal, 'k-');
title('input signal');
xlabel('Time (sec)');
ylabel('R*');

% Plot output signal at a particular (x, y) over time.
subplot(1,2,2);
outputSignal(1,:) = obj.ConeCurrentSignal(round(sz1/2),round(sz2/2),:);
plot((0:numel(outputSignal)-1)*dt, outputSignal, 'k-');
% if obj.noiseflag == 1
%     noisyOutputSignal = obj.ConeCurrentSignalPlusNoise(round(sz1/2),round(sz2/2),:);
%     plot((0:numel(outputSignal)-1)*dt, outputSignal, 'r-');
% end
title('output signal');
xlabel('Time (sec)');
ylabel('pA');
drawnow;