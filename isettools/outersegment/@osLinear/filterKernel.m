function filterKernel(obj,varargin)
% filterKernel: a method of @osLinear that generates the linear filters for
% the L, M and S cones that are used to convert isomerizations (R*) to
% outer segment current (pA). This method is called by the initialize
% method and the osLinearCompute method. The details of the generation of
% the impulse responses for the L, M and S cones are found in the
% filterConesLinear() utility function.
% 
    fprintf('<strong>\n%s:\n\t%s()\n</strong>', class(obj), mfilename());
    
    if length(varargin)==0
        [newIRFs, ~, ~] = filterConesLinear();
    else
        sensor = varargin{1};
        [newIRFs, ~, ~] = filterConesLinear(sensor);
    end
    
    obj = osLinearSet(obj, 'sconefilter', newIRFs(:,1),'units','pa');
    obj = osLinearSet(obj, 'mconefilter', newIRFs(:,2),'units','pa');
    obj = osLinearSet(obj, 'lconefilter', newIRFs(:,3),'units','pa');

end

