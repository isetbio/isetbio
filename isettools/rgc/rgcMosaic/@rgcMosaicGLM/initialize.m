function initialize(obj, varargin)
% intialize: a method of @rgcMosaicGLM that initializes the object
% following initialization by the superclass. This adds the generator
% function, the post spike filter and the coupling filters. This function
% is only called by rgcMosaicGLM, which itself is only called by rgcGLM or
% rgcCreate.
% 
%       rgcMosaicGLM.initialize(rgc, sensor, outersegment, varargin{:});
% 
% Inputs: 
%       rgc: an isetbio rgcGLM object
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
% See also rgcCreate, rgcGLM, rgcMosaicGLM.
% 
% (c) isetbio
% 09/2015 JRG

%     obj.generatorFunction = @erf;
    obj.generatorFunction = @exp;
%     obj.generatorFunction = @(x) 10*erf(x);
    
    obj.postSpikeFilter = buildPostSpikeFilter(.01);
    
    [obj.couplingFilter, obj.couplingMatrix] = buildCouplingFilters(obj, .01);



