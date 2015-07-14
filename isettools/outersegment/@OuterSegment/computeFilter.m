function computeFilter(obj)
    fprintf('<strong>\n%s:\n\t%s()\n</strong>', class(obj), mfilename());
    
%     obj.filterKernel = [];%exp(-obj.time/obj.timeConstant);

end

