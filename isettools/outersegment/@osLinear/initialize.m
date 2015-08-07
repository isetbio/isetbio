function initialize(obj)
% initialize: a method of @osLinear, initializes the osLinear object.
    fprintf('<strong>\n%s:\n\t%s()\n</strong>', class(obj), mfilename()); 
    filterKernel(obj);
%     sConeFilter = [];
%     mConeFilter = [];
%     lConeFilter = [];
end

