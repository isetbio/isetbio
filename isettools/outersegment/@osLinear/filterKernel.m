function filterKernel(obj,varargin)
% filterKernel: a method of @osLinear that generates the linear filters for
% the L, M and S cones that are used to convert isomerizations (R*) to
% outer segment current (pA). This method is called by the initialize
% method and the osLinearCompute method. The details of the generation of
% the impulse responses for the L, M and S cones are found in the
% filterConesLinear() utility function.
%
% Inputs: a cone osLinear object.
% Outputs: a cone osLinear object with the filter fields set to the temporal
% impulse response via the filterConesLinear function.
%
% 7/2015 JRG
% 
   
    if isempty(varargin)
        [newIRFs, ~, ~] = filterConesLinear();
    else
        sensor = varargin{1};
        [newIRFs, ~, ~] = filterConesLinear(sensor);
    end
    
    obj = osSet(obj, 'lconefilter', newIRFs(:,1),'units','pa');
    obj = osSet(obj, 'mconefilter', newIRFs(:,2),'units','pa');
    obj = osSet(obj, 'sconefilter', newIRFs(:,3),'units','pa');

end

