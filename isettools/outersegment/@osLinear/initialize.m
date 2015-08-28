function initialize(obj)
% initialize: a method of @osLinear, initializes the osLinear object.
% 
% Inputs: the osLinear object.
% Outputs: the osLinear object sConeFilter, mConeFilter and lConeFilter
% properties store the filter impulse response functions.
% 
% See the utility filterConeLinear for details of the generation of the
% impulse responses.
% 
% 7/2015 JRG


filterKernel(obj); % builds the L, M and S-cone linear temporal filters
    
end

