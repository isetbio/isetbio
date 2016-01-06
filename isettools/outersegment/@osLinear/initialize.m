function obj = initialize(obj,varargin)
% initialize: a method of @osLinear, initializes the osLinear object.
% 
% Inputs: the osLinear object.
% Outputs: the osLinear object sConeFilter, mConeFilter and lConeFilter
% properties store the filter impulse response functions.
% 
% See the utility osFilterConesLinear for details of the generation of the
% impulse responses.
% 
% 7/2015 JRG

if ~isempty(varargin)
    sensor = varargin{1};
    newIRFs = osFilterConesLinear(sensor);
else
    newIRFs = osFilterConesLinear();
end
obj = osSet(obj, 'lconefilter', newIRFs(:,1));
obj = osSet(obj, 'mconefilter', newIRFs(:,2));
obj = osSet(obj, 'sconefilter', newIRFs(:,3));

end

