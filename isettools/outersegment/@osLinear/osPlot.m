function h = osPlot(obj, absorptions, varargin)
% Plots about the outersegment properties
%
% Inputs: 
%   osLinear object
%   absorptions
% 
% Options:
%  isomerizations
%  filter kernels
%  current
%  all
%
% Outputs: 
%    h is a handle to the plot window
%
% Properties that can be plotted:
%
% Examples:
%   osL.plot(absorptions,'type','isomerizations')
%   osL.plot(absorptions,'type','current')
%   osL.plot(absorptions,'type','filter kernels')
%   osL.plot(absorptions,'type','all')
%
% (c) isetbio
% 09/2015 JRG

%% Check for the number of arguments and create parser object.

p = inputParser;
addRequired(p, 'obj');
addRequired(p, 'sensor');
addParameter(p,'type', 'all', @ischar);

p.parse(obj, absorptions, varargin{:});
params  = p.Results;
absorptions  = params.sensor;
type   = params.type;

% Choosing the plot type
switch ieParamFormat(type)
    
    case {'isomerizations'}
        % Isomerizations over time as input signals
        h = vcNewGraphWin;
        
        % Plot input signal (isomerizations) at a particular (x, y) over time.
        set(h, 'Name', sprintf('Isomerizations %s', class(obj)));       
        osPlotIsomerizations(absorptions);

        
    case {'filterkernels','filter','filters'}
        % Plot linear temporal filters for L, M and S cones.

        h = vcNewGraphWin;
        osPlotKernels(obj,absorptions);
        
    case {'current'}
        % Plot output signal at a particular (x, y) over time.

        h = vcNewGraphWin;
        osPlotCurrent(obj,absorptions);

    case{'all'}
        % Puts each of the main plots in a subplot within the window
        h = vcNewGraphWin([],'wide');
        
        subplot(1,3,1)
        osPlotIsomerizations(absorptions)
        
        subplot(1,3,2)
        osPlotKernels(obj,absorptions)
        
        subplot(1,3,3)
        osPlotCurrent(obj,absorptions)
        
    otherwise
        warning('Unknown plot type %s\n',type);
end

end


function osPlotIsomerizations(sensor)
% 
%
dt = sensorGet(sensor, 'time interval');

isomerizations1 = sensorGet(sensor,'photon rate');
[sz1, sz2, ~] = size(isomerizations1);
inputSignal = squeeze(isomerizations1(round(sz1/2),round(sz2/2),:));

plot((0:numel(inputSignal)-1)*dt, inputSignal, 'k-');
title('input signal');
xlabel('Time (sec)');
ylabel('R*/sec');

end

%
function osPlotKernels(obj,sensor)
%

dt = sensorGet(sensor, 'time interval');

hold on;
plot((0:numel(obj.lmsConeFilter(:, 3))-1)*dt, obj.lmsConeFilter(:, 3),'b');
plot((0:numel(obj.lmsConeFilter(:, 2))-1)*dt, obj.lmsConeFilter(:, 2),'g');
plot((0:numel(obj.lmsConeFilter(:, 1))-1)*dt, obj.lmsConeFilter(:, 1),'r');

title('L, M, S cone filter kernels');
xlabel('Time (sec)');
ylabel('pA / (R*/sec)');
hold off

end

%%
function osPlotCurrent(obj,sensor)
%

dt = sensorGet(sensor, 'time interval');

outputSignalTemp = osGet(obj,'cone current signal');
sz = osGet(obj,'size');
outputSignal = reshape(outputSignalTemp,sz(1)*sz(2),sz(3));

plot((0:size(outputSignal,2)-1)*dt, outputSignal(1+floor((sz(1)*sz(2)/100)*rand(200,1)),:));
title('Output current');
xlabel('Time (sec)');
ylabel('pA');

end

