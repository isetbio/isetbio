function [fullResponse, nlResponse] = fullConvolve(mosaic, spResponseCenter, spResponseSurround)
% fullConvolve: a util function of the @rgc parent class, for a separable
% STRF finds the 1D convolution of the temporal impulse response with the
% output signal of the spatial convolution operation.
% 
% Inputs:
% 
% Outputs:
% 
% Example:
% 
% (c) isetbio
% 09/2015 JRG

spResponseSize = size(spResponseCenter{1,1}(:,:,1,1));
nSamples = size(spResponseCenter{1,1},3);
channelSize = size(spResponseCenter{1,1},4);

nCells = size(mosaic.cellLocation);

fprintf('Temporal Convolution, %s:     \n', mosaic.cellType);

tic
for xcell = 1:nCells(1)
    for ycell = 1:nCells(2)
        
        for rgbIndex = 1:channelSize
            
            % fullResponseRSRGB = zeros(size(spResponseCenter{1,1}));
            
            temporalIRCenter = mosaic.tCenter{rgbIndex};
            temporalIRSurround = mosaic.tSurround{rgbIndex};
            
            spResponseCenterRS = reshape(squeeze(spResponseCenter{xcell,ycell}(:,:,:,rgbIndex)), spResponseSize(1)*spResponseSize(2), nSamples);
            spResponseSurroundRS = reshape(squeeze(spResponseSurround{xcell,ycell}(:,:,:,rgbIndex)), spResponseSize(1)*spResponseSize(2), nSamples);
    
            if (sum(temporalIRCenter-temporalIRSurround) == 0) 
                % if the temporal impulse responses for center and surround are the same, combine before convolution for efficiency                                             
                fullResponseRSCombined = convn(spResponseCenterRS-spResponseSurroundRS, temporalIRCenter','full');
                
                startPoint = length(temporalIRCenter)-1; endPoint = nSamples+length(temporalIRCenter)-1;
                fullResponseRSRGB(:,:,rgbIndex) = fullResponseRSCombined(:,startPoint:endPoint);
                                
            else            
                % if the temporal impulse responses for center and surround are different, convolve both               
                fullResponseRSCenter = convn(spResponseCenterRS, temporalIRCenter','full');
                fullResponseRSSurround = convn(spResponseSurroundRS, temporalIRSurround','full');
                
                startPoint = length(temporalIRCenter)-1; endPoint = nSamples+length(temporalIRCenter)-1;
                fullResponseRSRGB(:,:,rgbIndex) = fullResponseRSCenter(:,startPoint:endPoint) - fullResponseRSSurround(:,startPoint:endPoint);
                
            end            
        end      
                        
        if isa(mosaic,'rgcMosaicSubunit');
            
            genFunction = mosaicGet(mosaic, 'generatorFunction');
            fullResponseRS = sum(genFunction(fullResponseRSRGB),3);            
            fullResponse{xcell,ycell,1} = mean(fullResponseRS);            
            nlResponse{xcell,ycell} = (mean(fullResponseRS,1));
        else
            fullResponseRS = sum(fullResponseRSRGB,3);                     
            fullResponse{xcell,ycell,1} = mean(fullResponseRS);            
            genFunction = mosaicGet(mosaic, 'generatorFunction');
            nlResponse{xcell,ycell} = genFunction(mean(fullResponseRS,1));
        end
        
        % % % fullResponse for RGB
        % fullResponse{xcell,ycell,2} =  reshape(fullResponseRS, spResponseSize(1), spResponseSize(2), size(fullResponseRS,2));
                                      
    end
end
toc