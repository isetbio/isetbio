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

% Required for Pillow code
global RefreshRate
RefreshRate = 100;    
%% Loop on the mosaics in the inner retina
for ii = 1:length(ir.mosaic)
    
    switch class(ir.mosaic{ii})
        case {'rgcGLM','rgcSubunit'}
            % Call the Pillow code to generate spikes for the whole mosaic
            % using the coupled GLM
            
            % Modified
            % responseSpikes = computeSpikesGLM(ir.mosaic{ii,1});
            
            % Wrappers for adapting isetbio mosaic properties to Pillow code
            glminput = setGLMinput(ir.mosaic{ii}.responseLinear);
            glmprs = setGLMprs(ir.mosaic{ii});
            % Run Pillow code
            [responseSpikesVec, Vmem] = simGLMcpl(glmprs, glminput');
            cellCtr = 0;
            nCells = size(ir.mosaic{ii}.responseLinear);
            for xc = 1:nCells(1)
                for yc = 1:nCells(2)
                    cellCtr = cellCtr+1;
                    responseSpikes{yc,xc} = responseSpikesVec{1,cellCtr};
                end
            end

            % Set mosaic property
            ir.mosaic{ii} = mosaicSet(ir.mosaic{ii},'responseSpikes', responseSpikes);
            
        case 'rgcLNP'
            % Call another version of the Pillow code for spike response
            % responseSpikes = computeSpikesPSF(ir.mosaic{ii});
            
            % Wrappers for adapting isetbio mosaic properties to Pillow code
            glminput = setGLMinput(ir.mosaic{ii}.responseLinear);
            
            % Uses post spike filter
            glmprs = setPSFprs(ir.mosaic{ii});
            % No post spike filter - break into different subclass?
            % glmprs = setLNPprs(ir.mosaic{ii});
            
            % Run Pillow code
            [responseSpikesVec, Vmem] = simGLMcpl(glmprs, glminput');
            cellCtr = 0;
            nCells = size(ir.mosaic{ii}.responseLinear);
            for xc = 1:nCells(1)
                for yc = 1:nCells(2)
                    cellCtr = cellCtr+1;
                    responseSpikes{yc,xc} = responseSpikesVec{1,cellCtr};
                end
            end
            
            ir.mosaic{ii} = mosaicSet(ir.mosaic{ii},'responseSpikes', responseSpikes);
            
        otherwise
            error('The rgcMosaic object is a model without a spike response; choose LNP or GLM for spikes.');
    end
end


end



