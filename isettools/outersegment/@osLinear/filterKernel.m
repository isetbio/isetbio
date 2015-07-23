function filterKernel(obj)
    fprintf('<strong>\n%s:\n\t%s()\n</strong>', class(obj), mfilename());
    
    [newIRFs, ~, ~] = filterConesLinear();
    obj.sConeFilter = newIRFs(:,1);
    obj.mConeFilter = newIRFs(:,2);
    obj.lConeFilter = newIRFs(:,3);

end

