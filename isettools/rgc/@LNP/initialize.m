function initialize(obj)
% intialize: a method of @rgcMosaicLNP that initializes the object
% following initialization by the superclass. This adds the generator
% function, the post spike filter and the coupling filters. This function
% is only called by rgcMosaicLNP, which itself is only called by rgcLNP or
% rgcCreate.
% 
%       rgcMosaicLNP.initialize(rgc, sensor, outersegment, varargin{:});
% 
% Inputs: 
%       rgc: an isetbio rgcLNP object
%       scene: an isetbio scene structure
%       sensor: an isetbio sensor structure
%       os: an isetbio outer segment structure
% 
% Outputs: the mosaic object with the generatorFunction, postSpikeFilter and
% couplingFilter properties set to appropriate values.
% 
% Example:
% 
%       obj.initialize(rgc, sensor, outersegment, varargin{:});
% 
% See also rgcCreate, rgcLNP, rgcMosaicLNP.
% 
% (c) isetbio
% 09/2015 JRG


% %% Add generator function
% % Need to make this into Gaussian CDF
% for cellTypeInd = 1:obj.numberCellTypes
    obj.generatorFunction = @exp;
    obj.numberTrials = 10;
% end

