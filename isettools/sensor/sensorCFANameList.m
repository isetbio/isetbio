function cfaName = sensorCFANameList
% List the valid ISET names for color filter arrays
%
%    cfaName = sensorCFAName(ISA)
%
% At this time, we support a short list of color filter array types.
% This returns the valid list, and it is used to define the popup boxes in
% sensorImageWindow. The list will grow in the future.
%   
% Examples:
%   cfaNames = sensorCFAName
%
% Copyright ImagEval Consultants, LLC, 2005

cfaName = ...
    {'bayer-grbg','bayer-rggb','bayer-bggr','bayer-gbrg',...
    'bayer-ycmy','Monochrome','Four Color','Other'}; 

return; 

