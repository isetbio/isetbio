function obj = rgcCompute(obj, outerSegment, varargin)
% rgcCompute: a method of @rgc that computes the spiking output of the
% rgc mosaic to an arbitrary stimulus.
% 


obj = rgcCompute@rgc(obj, outerSegment, varargin{:});

% Linear and nonlinear responses calculated in @rgc/rgcCompute

for cellTypeInd = 1:length(obj.mosaic)
            
    spikeResponse = computeSpikesGLM(obj.mosaic{cellTypeInd,1});   
    obj.mosaic{cellTypeInd} = mosaicSet(obj.mosaic{cellTypeInd},'spikeResponse', spikeResponse);
    
    [raster, psth] = computePSTH(obj.mosaic{cellTypeInd,1});
    
    obj.mosaic{cellTypeInd} = mosaicSet(obj.mosaic{cellTypeInd},'rasterResponse',raster);
    obj.mosaic{cellTypeInd} = mosaicSet(obj.mosaic{cellTypeInd},'psthResponse',psth);
        
    clear spikeResponse raster psth
end




