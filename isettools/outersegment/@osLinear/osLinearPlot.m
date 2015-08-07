% Method to plot the computed results of the outsergment object
function plot(obj, sensor)

dt = sensorGet(sensor, 'time interval');

fprintf('<strong>\n%s:\n\t%s()\n</strong>', class(obj), mfilename());
figNum = 1; 
h = figure();
set(h, 'Name', sprintf('Output of %s', class(obj)));
set(h, 'Position', [10+50*figNum 10+50*figNum, 1024 256]);
subplot(1,3,1);
[sz1 sz2 sz3] = size(sensor.data.volts); 
inputSignal(1,:) = sensor.data.volts(round(sz1/2),round(sz2/2),:);
plot((0:numel(inputSignal)-1)*dt, inputSignal, 'k-');
title('input signal');
xlabel('Time (sec)');
ylabel('R*');
subplot(1,3,2);
hold on;
plot((0:numel(obj.sConeFilter)-1)*dt, obj.sConeFilter,'b'); 
plot((0:numel(obj.mConeFilter)-1)*dt,obj.mConeFilter,'g'); 
plot((0:numel(obj.lConeFilter)-1)*dt,obj.lConeFilter,'r');
title('L, M, S cone filter kernels');
xlabel('Time (sec)');
ylabel('pA');
subplot(1,3,3);
outputSignal(1,:) = obj.ConeCurrentSignal(round(sz1/2),round(sz2/2),:);
plot((0:numel(outputSignal)-1)*dt, outputSignal, 'k-');
title('output signal');
xlabel('Time (sec)');
ylabel('pA');
drawnow;