function ir = irComputeSpikes(ir, varargin)
% Generate spikes for certain rgc objects
%
% We have models that convert from the continuous response to the spiking
% output for several types of rgc mosaics.  This is the gateway routine
% that examines which rgc type we have and invokes the proper continuous to
% spike response for that type of rgc.
%
% Inputs: the inner retina object
%
% Outputs: the spikes responses are attached to the rgc mosaics in the ir
% object.
%
% Example:
%  os = osCreate('identity');
%  ir = irCreate(os,'model','glm');
%  ir.mosaicCreate('mosaicType','on midget');
%  ir.computeContinuous;
%  ir.computeSpike;
%
% JRG (c) isetbio


%% Loop on the mosaics in the inner retina
for ii = 1:length(ir.mosaic)
    
    switch class(ir.mosaic{ii})
        case {'rgcGLM','rgcSubunit'}
            % Call the Pillow code to generate spikes for the whole mosaic
            % using the coupled GLM
            responseSpikes = computeSpikesGLM(ir.mosaic{ii,1});
            ir.mosaic{ii} = mosaicSet(ir.mosaic{ii},'responseSpikes', responseSpikes);
            
        case 'rgcLNP'
            % Call another version of the Pillow code for spike response
            responseSpikes = computeSpikesPSF(ir.mosaic{ii});
            ir.mosaic{ii} = mosaicSet(ir.mosaic{ii},'responseSpikes', responseSpikes);
            
        otherwise
            error('The rgcMosaic object is a model without a spike response; choose LNP or GLM for spikes.');
    end
end


end



