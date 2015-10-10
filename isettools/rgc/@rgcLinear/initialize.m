function initialize(obj, sensor, outersegment, varargin)
% intialize: a method of @rgcLNP that initializes the object based on a
% series of input parameters that can include the location of the
% retinal patch.
% 
% Inputs:
% 
% Outputs:
% 
% Example:
% 
% (c) isetbio
% 09/2015 JRG


for cellTypeInd = 1:length(obj.mosaic)
    obj.mosaic{cellTypeInd} = rgcMosaicLinear(cellTypeInd, obj, sensor, outersegment, varargin{:});
end