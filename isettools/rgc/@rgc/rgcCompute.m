function obj = rgcCompute(obj, outerSegment, varargin)
% rgcCompute: a method of @rgc that computes the spiking output of the
% rgc mosaic to an arbitrary stimulus.
% 
% Inputs:
% 
% Outputs:
% 
% Example:
% 
% (c) isetbio
% 09/2015 JRG

%% 

for cellTypeInd = 3%1:obj.numberCellTypes
    
    nCells = size(obj.mosaic{cellTypeInd}.spatialRFArray);
    
    if isa(outerSegment,'osIdentity')
        spConvolveStimulus = osGet(outerSegment, 'rgbData');
%         spConvolveStimulus = spConvolveStimulus - mean(spConvolveStimulus(:));
%         spConvolveStimulus = 10*spConvolveStimulus./max(abs(spConvolveStimulus(:)));
    else
        spConvolveStimulus = osGet(outerSegment, 'coneCurrentSignal');
        spConvolveStimulus = spConvolveStimulus - mean(spConvolveStimulus(:));
        spConvolveStimulus = 5*spConvolveStimulus./max(spConvolveStimulus(:));
    end   
    
    % Given separable STRF, convolve 2D spatial RF then 1D temporal response.
     
    spResponse = spConvolve(obj.mosaic{cellTypeInd,1}, spConvolveStimulus);
       
%     if isa(outerSegment,'osIdentity');    
        [fullResponse, nlResponse] = fullConvolve(obj.mosaic{cellTypeInd,1}, spResponse);        
%     else        
%         
%         genFunc = obj.mosaic{cellTypeInd}.generatorFunction;
%         for xcell = 1:nCells(1)
%             for ycell = 1:nCells(2)
%                 fullResponse{xcell,ycell} = mean((spResponse{xcell,ycell,1}) - (spResponse{xcell,ycell,2}),3);
%                 nlResponse{xcell,ycell} = genFunc(fullResponse{xcell,ycell});
%             end
%         end
%     end
        
    
    obj.mosaic{cellTypeInd} = mosaicSet(obj.mosaic{cellTypeInd},'linearResponse', fullResponse);
    
    if ~isa(obj, 'rgcLinear')
        
        obj.mosaic{cellTypeInd} = mosaicSet(obj.mosaic{cellTypeInd},'nlResponse', nlResponse);
        
        %%% Compute spikes
        
        % for cellTypeInd = 1%:4%obj.numberCellTypes
        %     obj.mosaic{cellTypeInd,1}.spikeResponse = computeSpikes(obj.mosaic{cellTypeInd,1}.nlResponse, sensor, outersegment);
        fprintf('Spike generation, %s:      \n', obj.mosaic{cellTypeInd}.nameCellType);
        tic
        
        if isa(obj,'rgcLNP');
            spikeResponse = computeSpikes(obj.mosaic{cellTypeInd,1}.nlResponse, obj.mosaic{cellTypeInd}.postSpikeFilter);
            % elseif 0
            %     spikeResponse = computeSpikesPSF(obj.mosaic{cellTypeInd,1}.nlResponse, obj.mosaic{cellTypeInd}.postSpikeFilter, sensor, outersegment);
        elseif isa(obj,'rgcGLM')
            spikeResponse = computeSpikesGLM(obj.mosaic{cellTypeInd,1});
        end
        % obj = rgcMosaicSet(obj, 'spikeResponse', spikeResponse);
        
        obj.mosaic{cellTypeInd} = mosaicSet(obj.mosaic{cellTypeInd},'spikeResponse', spikeResponse);
        
        [raster psth] = computePSTH(obj.mosaic{cellTypeInd,1});
        
        obj.mosaic{cellTypeInd} = mosaicSet(obj.mosaic{cellTypeInd},'rasterResponse',raster);
        obj.mosaic{cellTypeInd} = mosaicSet(obj.mosaic{cellTypeInd},'psthResponse',psth);
        
        toc
        
    end
    ph=1;
end
close;
