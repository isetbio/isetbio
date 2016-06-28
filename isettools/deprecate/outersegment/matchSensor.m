function obj = matchSensor(obj,varargin)
% Matches the size of the osLinear object to the sensor data
%
% Inputs: the osLinear object.
%
% Outputs: the osLinear object sConeFilter, mConeFilter and lConeFilter
% properties store the filter impulse response functions.
% 
% See the utility osFilterConesLinear for details of the generation of the
% impulse responses.
% 
% 7/2015 JRG

if ~isempty(varargin)
    % There is a sensor.  Match the impulse response functions
    sensor = varargin{1};
else
    % Uses the defaults
    sensor = [];
end
newIRFs = osFilterConesLinear(sensor);

% Add the IRFs to the object
obj = osSet(obj, 'lconefilter', newIRFs(:,1));
obj = osSet(obj, 'mconefilter', newIRFs(:,2));
obj = osSet(obj, 'sconefilter', newIRFs(:,3));

end

