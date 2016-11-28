function osPlot(obj, absorptions, varargin)
% Plots the input (photons/sec), linear filters and output (pA) of the
% linear outer segment.
%
% Inputs: 
%   osLinear object
%   absorptions
% 
% Options:
%  isomerizations
%  current
%  all
%
% Outputs: 
%    h is a handle to the plot window
%
% Properties that can be plotted:
%
% Examples:
%   osBP.plot(absorptions,'type','isomerizations')
%   osBP.plot(absorptions,'type','current')
%   osBP.plot(absorptions,'type','all')
%
% (c) isetbio
% 09/2015 JRG

%% Check for the number of arguments and create parser object.

p = inputParser;
addRequired(p, 'obj');
addRequired(p, 'sensor');
addParameter(p,'type', 'all', @ischar);

p.parse(obj, absorptions, varargin{:});
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
        
    case {'current'}
        % Plot output signal at a particular (x, y) over time.
        h = vcNewGraphWin;
        osPlotCurrent(obj,absorptions);

    case{'all'}
        % Puts each of the main plots in a subplot within the window
        h = vcNewGraphWin([],'wide');
        
        subplot(1,2,1)
        osPlotIsomerizations(absorptions)
       
        subplot(1,2,2)
        osPlotCurrent(obj,absorptions)
        
    otherwise
        % Pass back to super class plot routine
        plot@outersegment(obj,param,varargin{:});
        
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

