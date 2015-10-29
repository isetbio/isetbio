function initialize(obj, varargin)
% intialize: a method of @rgcMosaicGLM that initializes the object
% following initialization by the superclass. This adds the generator
% function, the post spike filter and the coupling filters.
% 
% Inputs: the rgcGLM object.
% 
% Outputs: the object with the generatorFunction, postSpikeFilter and
% couplingFilter properties set to appropriate values.
% 
% Example:
% rgc1 = rgcCreate('glm', scene, sensor, os, 'right', 3.0, 180);
% 
% (c) isetbio
% 09/2015 JRG

%     obj.generatorFunction = @erf;
    obj.generatorFunction = @exp;
%     obj.generatorFunction = @(x) 10*erf(x);
    
    obj.postSpikeFilter = buildPostSpikeFilter(.01);
    
    [obj.couplingFilter, obj.couplingMatrix] = buildCouplingFilters(obj, .01);



