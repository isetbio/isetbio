function h = osPlot(obj, pType, varargin)
% Plots for the os object 
%
% Inputs: 
%   absorptions - input absorptions
%   filters     - impulse response filters
%   output
%   all
% 
% Outputs: plot(s)
% 
% Properties that can be plotted:
% 
% Examples:
%   osPlot(os, sensor);
%   osPlot(os, sensor,'input');
%   osPlot(os, sensor,'filters');
% 
% (c) isetbio
% 09/2015 JRG

%% Check for the number of arguments and create parser object.
p = inputParser; 
p.CaseSensitive = false; 
p.FunctionName = mfilename;

% Make key properties that can be set required arguments, and require
% values along with key names.
allowableFieldsToSet = {...
        'input'...
        'filter',...
        'filters',...
        'output',...
        'cell',...
        'all'...
    };

% vFunc = @(x)(isstruct(x) && ....)
p.addParameter('absorptions',[],@isstruct)
p.addOptional('value',@isnumeric);
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
    case{'input','absorptions'}
    
        
        % Plot input signal (isomerizations) at a particular (x, y) over time.
        h = vcNewGraphWin;
        dt = obj.get('time interval');
        
        % since data is in (x, y, t) format, choose an (x, y) value to observe over
        % timesubplot(1,3,1);
        sz = os.get('array size');
        isomerizations = os.get('isomerizations');
        isomerizations = squeeze(isomerizations(round(sz(1)/2),round(sz(2)/2),:));
        
        plot((0:numel(isomerizations)-1)*dt, isomerizations, 'k-');
        title('input signal');
        xlabel('Time (sec)');
        ylabel('R*');
        
    case{'filter','filters'}
        
        dt = sensorGet(sensor, 'time interval');
        
        % Plot input signal (isomerizations) at a particular (x, y) over time.
        vcNewGraphWin;
        
        % Plot linear temporal filters for L, M and S cones.
        
        hold on;
        plot((0:numel(obj.sConeFilter)-1)*dt, obj.sConeFilter,'b');
        plot((0:numel(obj.mConeFilter)-1)*dt, obj.mConeFilter,'g');
        plot((0:numel(obj.lConeFilter)-1)*dt, obj.lConeFilter,'r');
        title('L, M, S cone filter kernels');
        xlabel('Time (sec)');
        ylabel('pA / (R*/sec)');
        
    case{'output'}
        % Cone output current
        h = vcNewGraphWin;
        
        dt = obj.get('time interval');

        % Plot current signal at a particular (x, y) over time
        sz = os.get('array size');
        outputSignal(1,:) = obj.ConeCurrentSignal(round(sz(1)/2),round(sz(2)/2),:);
        plot((0:numel(outputSignal)-1)*dt, outputSignal, 'k-');
        title('output signal');
        xlabel('Time (sec)');
        ylabel('pA');
        
    case{'all'}
        % Recommend writing this as three calls to the stuff above
        
        dt = sensorGet(sensor, 'time interval');
                
        % Plot input signal (isomerizations) at a particular (x, y) over time.
        h = vcNewGraphWin([],'wide');
        set(h, 'Name', sprintf('Output of %s', class(obj)));
        
        % since data is in (x, y, t) format, choose an (x, y) value to observe over
        % timesubplot(1,3,1);
        subplot(1,3,1)
        isomerizations1 = sensorGet(sensor,'photons');
        [sz1, sz2, sz3] = size(isomerizations1);
        inputSignal = squeeze(isomerizations1(round(sz1/2),round(sz2/2),:));
        plot((0:numel(inputSignal)-1)*dt, inputSignal, 'k-');
        title('input signal');
        xlabel('Time (sec)');
        ylabel('R*');
        
        % Plot linear temporal filters for L, M and S cones.
        subplot(1,3,2);
        hold on;
        plot((0:numel(obj.sConeFilter)-1)*dt, obj.sConeFilter,'b');
        plot((0:numel(obj.mConeFilter)-1)*dt, obj.mConeFilter,'g');
        plot((0:numel(obj.lConeFilter)-1)*dt, obj.lConeFilter,'r');
        title('L, M, S cone filter kernels');
        xlabel('Time (sec)');        
        ylabel('pA / (R*/sec)');
        
        % Plot output signal at a particular (x, y) over time.
        subplot(1,3,3);
        outputSignalTemp = osGet(obj,'coneCurrentSignal');
        if isfield(params, 'cell');
            outputSignal(1,:) = outputSignalTemp(params.cell(1),params.cell(2),:);
        else
            % outputSignal(1,:) = outputSignalTemp(round(sz1/2),round(sz2/2),:);
            outputSignal = reshape(outputSignalTemp,sz1*sz2,sz3);
        end
        plot((0:size(outputSignal,2)-1)*dt, outputSignal(1+floor((sz1*sz2/100)*rand(200,1)),:));
        title('output signal');
        xlabel('Time (sec)');
        ylabel('pA');

        
end
