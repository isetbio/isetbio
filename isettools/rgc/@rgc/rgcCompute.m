function obj = rgcCompute(obj, sensor, outersegment, varargin)
% rgcCompute: a method of @rgcLNP that computes the spiking output of the
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

%% Linear output via separable convolution 

for cellTypeInd = 1:obj.numberCellTypes
    % spConvolveFilter = rgcGet(obj, 'mosaic', cellTypeInd);
%     spConvolveStimulus = osGet(outersegment, 'coneCurrentSignal');
    spConvolveStimulus3Chann = osGet(outersegment, 'rgbData');
    
    % NEED TO FIX THIS THING!!!
    % ONLY DOING ONE CHANNEL OUT OF RGB
    % spConvolveStimulus = squeeze(spConvolveStimulus3Chann(:,:,1,:));
    spConvolveStimulus = spConvolveStimulus3Chann;
    
    % Given separable STRF, first convolve 2D spatial RF then 1D temporal
    % response.
    %  spResponse{cellTypeInd,1} = spConvolve(spConvolveFilter, spConvolveStimulus);
    spResponse = spConvolve(obj.mosaic{cellTypeInd,1}, spConvolveStimulus);
   
    
    [fullResponse, nlResponse] = fullConvolve(obj.mosaic{cellTypeInd,1}, spResponse);
    
    obj.mosaic{cellTypeInd} = mosaicSet(obj.mosaic{cellTypeInd},'linearResponse', fullResponse);
%     obj.mosaic{1}.set('linearResponse', fullResponse);

if ~isa(obj, 'rgcLinear')
    
    obj.mosaic{cellTypeInd} = mosaicSet(obj.mosaic{cellTypeInd},'nlResponse', nlResponse);
    
   
    
%     obj.mosaic{cellTypeInd,1}.linearResponse = fullResponse;
%     obj.mosaic{cellTypeInd,1}.nlResponse = nlResponse;

%     toc
%     tic
%     % Apply generator function
%     obj.mosaic{cellTypeInd,1}.nlResponse = generatorLookup(fullResponse);
%     toc
% end

%% Generator function
% 
% for cellTypeInd = 1%:obj.numberCellTypes        
%     obj.mosaic{cellTypeInd,1}.nlResponse = generatorLookup(obj.mosaic{cellTypeInd,1}.linearResponse);    
% end

%% Compute spikes

% for cellTypeInd = 1%:4%obj.numberCellTypes        
%     obj.mosaic{cellTypeInd,1}.spikeResponse = computeSpikes(obj.mosaic{cellTypeInd,1}.nlResponse, sensor, outersegment); 
    fprintf('Spike generation, %s:      \n', obj.mosaic{cellTypeInd}.nameCellType);
    tic
    
    if 0    
        spikeResponse = computeSpikes(obj.mosaic{cellTypeInd,1}.nlResponse, obj.mosaic{cellTypeInd}.postSpikeFilter, sensor, outersegment);   
    elseif 0
        spikeResponse = computeSpikesPSF(obj.mosaic{cellTypeInd,1}.nlResponse, obj.mosaic{cellTypeInd}.postSpikeFilter, sensor, outersegment); 
    elseif 1       
        spikeResponse = computeSpikesGLM(obj.mosaic{cellTypeInd,1}, sensor, outersegment);       
    end
%     obj = rgcMosaicSet(obj, 'spikeResponse', spikeResponse);
    
    obj.mosaic{cellTypeInd} = mosaicSet(obj.mosaic{cellTypeInd},'spikeResponse', spikeResponse);
    toc
    
end
end
