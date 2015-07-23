% Initialize LinearOuterSegment object
function initialize(obj)
    fprintf('<strong>\n%s:\n\t%s()\n</strong>', class(obj), mfilename()); 
    filterKernel(obj);
%     sConeFilter = [];
%     mConeFilter = [];
%     lConeFilter = [];
end

