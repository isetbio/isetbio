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
if nargin == 2; varargin{1} = 'output'; end;
p = inputParser; p.CaseSensitive = false; p.FunctionName = mfilename;

% This flag causes the parser not to throw an error here in the superclass
% call. The subclass call will throw an error.
% p.KeepUnmatched = true;

% Make key properties that can be set required arguments, and require
% values along with key names.
allowableFieldsToSet = {...
    'output'...
       
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

    case{'output'}
        
        dt = sensorGet(sensor, 'time interval');
        
        % Plot input signal (isomerizations) at a particular (x, y) over time.
        h = vcNewGraphWin;
        
        % Plot output signal at a particular (x, y) over time.
        [sz1 sz2 sz3] = size(obj.photonRate);
        outputSignal = squeeze(obj.photonRate(round(sz1/2),round(sz2/2),:,:));
        plot((0:numel(outputSignal(:,1))-1)*dt, outputSignal, 'k-');
        
        title('output signal, RGB');
        xlabel('Time (sec)');
        ylabel('Luminance');
        
end
