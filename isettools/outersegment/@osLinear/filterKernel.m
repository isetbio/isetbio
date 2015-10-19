function obj = filterKernel(obj,varargin)
% Generates the L,M, and S linear temporal filters to convert isomerizations (R*) to
% outer segment current (pA).
%
% This method is called by the initialize method and the osLinearCompute
% method. The details of the generation of the impulse responses for the L,
% M and S cones are found in the filterConesLinear() utility function.
%
% The time base of the filters is set from the sensor object ...
% Let's discuss how that should be denoted here.
%
% Inputs: a cone osLinear object.
% Outputs: a cone osLinear object with the filter fields set to the temporal
% impulse response via the filterConesLinear function.
%
% 7/2015 JRG

sensor = [];
if ~isempty(varargin), sensor = varargin{1}; end

newIRFs = filterConesLinear(sensor);

obj = osSet(obj, 'lconefilter', newIRFs(:,1),'units','pa');
obj = osSet(obj, 'mconefilter', newIRFs(:,2),'units','pa');
obj = osSet(obj, 'sconefilter', newIRFs(:,3),'units','pa');

end

