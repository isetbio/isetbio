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
    
%     obj = obj.osLinearSet('sconefilter', newIRFs(:,1));
%     obj = obj.osLinearSet('mConeFilter', newIRFs(:,2));
%     obj = obj.osLinearSet('lConeFilter', newIRFs(:,3));

    obj = osLinearSet(obj, 'sconefilter', newIRFs(:,1),'units','pa');
    obj = osLinearSet(obj, 'mconefilter', newIRFs(:,2),'units','pa');
    obj = osLinearSet(obj, 'lconefilter', newIRFs(:,3),'units','pa');

end

