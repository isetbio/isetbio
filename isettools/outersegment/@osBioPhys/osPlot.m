function osPlot(obj, sensor, varargin)
% osPlot: a method of @oueterSegment that plots os object 
% properties using the input parser structure.
% 
% Inputs: os object, sensor, property to be plotted
% 
% Outputs: plot(s)
% 
% Properties that can be plotted:
% 
% Examples:
%   osPlot(os, sensor);
%   osPlot(os, sensor,'input');
%   osPlot(os, sensor,'output');
% 
% (c) isetbio
% 09/2015 JRG

% Check for the number of arguments and create parser object.
% Parse key-value pairs.
% 
% Check key names with a case-insensitive string, errors in this code are
% attributed to this function and not the parser object.
error(nargchk(0, Inf, nargin));
% if there is no argument for the type of plot, set default to all:
if nargin == 2; varargin{1} = 'all'; end;
p = inputParser; p.CaseSensitive = false; p.FunctionName = mfilename;

% This flag causes the parser not to throw an error here in the superclass
% call. The subclass call will throw an error.
% p.KeepUnmatched = true;

% Make key properties that can be set required arguments, and require
% values along with key names.
allowableFieldsToSet = {...
        'input'...
        'output',...
        'all'...
    };
p.addRequired('what',@(x) any(validatestring(x,allowableFieldsToSet)));

% % Define what units are allowable.
% allowableUnitStrings = {'a', 'ma', 'ua', 'na', 'pa'}; % amps to picoamps
% 
% % Set up key value pairs.
% % Defaults units:
% p.addParameter('units','pa',@(x) any(validatestring(x,allowableUnitStrings)));

% Parse and put results into structure p.
p.parse(varargin{:}); params = p.Results;

% Set key-value pairs.
switch lower(params.what)
    case{'input'}
    
        dt = sensorGet(sensor, 'time interval');
        
        % Plot input signal (isomerizations) at a particular (x, y) over time.
        h = vcNewGraphWin;
        
        % since data is in (x, y, t) format, choose an (x, y) value to observe over
        % timesubplot(1,3,1);
        
        isomerizations1 = sensorGet(sensor,'photons');
        [sz1 sz2 sz3] = size(isomerizations1);
        inputSignal = squeeze(isomerizations1(round(sz1/2),round(sz2/2),:));
        plot((0:numel(inputSignal)-1)*dt, inputSignal, 'k-');
        title('input signal');
        xlabel('Time (sec)');
        ylabel('R*');

    case{'output'}
        % Need to allow passing in which pixel or even an ROI
        
        dt = sensorGet(sensor, 'time interval');
        
        % Plot input signal (isomerizations) at a particular (x, y) over time.
        h = vcNewGraphWin;
        
        sz= sensorGet(sensor,'size');
        
        % Plot output signal at a particular (x, y) over time.
        
        outputSignal(1,:) = obj.ConeCurrentSignal(round(sz(1)/2),round(sz(2)/2),:);
        plot((0:numel(outputSignal)-1)*dt, outputSignal, 'k-');
        title('output signal');
        xlabel('Time (sec)');
        ylabel('pA');
        
    case{'all'}
        
        dt = sensorGet(sensor, 'time interval');
                
        % Plot input signal (isomerizations) at a particular (x, y) over time.
        h = vcNewGraphWin([],'wide');
        set(h, 'Name', sprintf('Output of %s', class(obj)));
        
        
        % Plot input signal (isomerizations) at a particular (x, y) over time.
        subplot(1,2,1);
        isomerizations1 = sensorGet(sensor,'photons');
        [sz1 sz2 sz3] = size(isomerizations1);
        inputSignal = squeeze(isomerizations1(round(sz1/2),round(sz2/2),:));
        plot((0:numel(inputSignal)-1)*dt, inputSignal, 'k-');
        title('input signal');
        xlabel('Time (sec)');
        ylabel('R*');
        
        % Plot output signal at a particular (x, y) over time.
        subplot(1,2,2);        
        outputSignalTemp = osGet(obj,'coneCurrentSignal');
        outputSignal(1,:) = outputSignalTemp(round(sz1/2),round(sz2/2),:);
        plot((0:numel(outputSignal)-1)*dt, outputSignal, 'k-');
        title('output signal');
        xlabel('Time (sec)');
        ylabel('pA');
        
end
