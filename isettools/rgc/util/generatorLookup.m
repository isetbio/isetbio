function nlResponse = generatorLookup(linearResponse)
% generatorFunction: a util function of the @rgc parent class, this
% converts the linear response of the STRF to a nonlinear response that
% maps more closely onto real RGC firing. The generator function is a
% lookup table that maps the linear response onto the neuron's continuous
% firing rate.
% 
% Inputs:
% 
% Outputs:
% 
% Example:
% 
% (c) isetbio
% 09/2015 JRG


spResponseSize = size(linearResponse{1,1}(:,:,1));
nSamples = size(linearResponse{1,1},3);

nCells = size(linearResponse);

for xcell = 1:nCells(1)
    for ycell = 1:nCells(2)
        % NEED TO ALLOW USER TO SPECIFY GENERATOR, ALSO CDF OF GAUSSIAN
        nlResponse{xcell,ycell} = exp(mean(mean(linearResponse{xcell,ycell},1),2));
        
    end
end

