function osPlot(obj, sensor, varargin)
% Plots the input (photons/sec), linear filters and output (pA) of the
% linear outer segment.
%
% Inputs: os object, sensor, property to be plotted
%
% Outputs: plot(s)
%
% Properties that can be plotted:
%
% Examples:
%   osPlot(os, sensor);
%
% (c) isetbio
% 09/2015 JRG

%% Check for the number of arguments and create parser object.

p = inputParser;
addRequired(p, 'obj');
addRequired(p, 'sensor');
% addParameter(p, 'sensor', 'sensor', @isstruct);
addParameter(p, 'type', 'all', @isstring);

p.parse(obj, sensor, varargin{:});

params = p.Results;
sensor = params.sensor;
type = params.type;

% Set key-value pairs.
switch ieParamFormat(params.type)
    
    case{'all'}
        if exist('sensor','var')
            dt = sensorGet(sensor, 'time interval');
            
            
            % Plot input signal (isomerizations) at a particular (x, y) over time.
            h = vcNewGraphWin([],'wide');
            set(h, 'Name', sprintf('Output of %s', class(obj)));
            
            % since data is in (x, y, t) format, choose an (x, y) value to observe over
            % timesubplot(1,3,1);
            subplot(1,3,1)
            isomerizations1 = sensorGet(sensor,'photons');
            [sz1 sz2 sz3] = size(isomerizations1);
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
        else
            warning('Need sensor input');
        end
        
    otherwise
        warning('Need sensor input');
end
