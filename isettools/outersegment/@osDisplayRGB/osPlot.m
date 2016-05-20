function osPlot(obj, sensor, varargin)
% Plots the input (photons/sec), linear filters and output (pA) of the
% linear outer segment.
%
% Inputs: os object, sensor
%
% Outputs: plot(s)
%
% Properties that can be plotted:
%         Output
%
% Examples:
%   osPlot(os, sensor);
%
% (c) isetbio
% 09/2015 JRG

%%
p = inputParser;
addRequired(p, 'obj');
addRequired(p, 'sensor');
% addParameter(p, 'sensor', 'sensor', @isstruct);
addParameter(p, 'type', 'all', @isstring);

p.parse(obj, sensor, varargin{:});

params = p.Results;
sensor = params.sensor;

% Set key-value pairs.
switch ieParamFormat(params.type)

    case{'output'}
        
        dt = sensorGet(sensor, 'time interval');
        
        % Plot input signal (isomerizations) at a particular (x, y) over time.
        h = vcNewGraphWin;
        
        % Plot output signal at a particular (x, y) over time.
        [sz1 sz2 sz3] = size(obj.rgbData);
        outputSignal = squeeze(obj.rgbData(round(sz1/2),round(sz2/2),:,:));
        plot((0:numel(outputSignal(:,1))-1)*dt, outputSignal, 'k-');
        
        title('output signal, RGB');
        xlabel('Time (sec)');
        ylabel('Luminance');
        
end
