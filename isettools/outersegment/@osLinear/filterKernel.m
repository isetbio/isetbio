function filterKernel(obj,varargin)
    fprintf('<strong>\n%s:\n\t%s()\n</strong>', class(obj), mfilename());
    
    if length(varargin)==0
        [newIRFs, ~, ~] = filterConesLinear();
    else
        sensor = varargin{1};
        [newIRFs, ~, ~] = filterConesLinear(sensor);
    end
    
%     obj.sConeFilter = newIRFs(:,1);
%     obj.mConeFilter = newIRFs(:,2);
%     obj.lConeFilter = newIRFs(:,3);
    
    obj = obj.osLinearSet('sConeFilter', newIRFs(:,1));
    obj = obj.osLinearSet('mConeFilter', newIRFs(:,2));
    obj = obj.osLinearSet('lConeFilter', newIRFs(:,3));

end

