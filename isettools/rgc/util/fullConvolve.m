function [fullResponse, nlResponse] = fullConvolve(mosaic, spResponse)
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

spResponseSize = size(spResponse{1,1}(:,:,1));
nSamples = size(spResponse{1,1},4);

nCells = size(mosaic.spatialRFArray);

rfSize = size(mosaic.spatialRFArray{1,1});

% fullResponse = cell(nCells); nlResponse = cell(nCells);
fprintf('Temporal Convolution, %s:     \n', mosaic.nameCellType);


%         tempIRmatrix = convmtx(temporalIR,size(spResponseRS,2));


%     fprintf('RGB = %d     \n', rgbIndex);
    tic
for xcell = 1:nCells(1)
    for ycell = 1:nCells(2) 
        
for rgbIndex = 1:3


            if 0%mosaic.temporalImpulseResponseCenterRGB == mosaic.temporalImpulseResponseSurroundRGB
                spResponseRS = reshape(squeeze(spResponse{xcell,ycell}(:,:,rgbIndex,:)), spResponseSize(1)*spResponseSize(2), nSamples);
            
                temporalIR = mosaic.temporalImpulseResponseCenterRGB{rgbIndex};
                
%                  tic
%                 tempIRmatrix = convmtx(temporalIR,size(spResponseRS,2));
%                 fullResponseRS = (tempIRmatrix*spResponseRS');
%                 toc
                
%                 tic
                fullResponseRSRGB(:,:,rgbIndex) = convn(spResponseRS, temporalIR');
%                 toc
                
               
            elseif 1
                
                %         tic
                
                spResponseCenterRS = reshape(squeeze(spResponse{xcell,ycell,1}(:,:,rgbIndex,:)), spResponseSize(1)*spResponseSize(2), nSamples);
                spResponseSurroundRS = reshape(squeeze(spResponse{xcell,ycell,2}(:,:,rgbIndex,:)), spResponseSize(1)*spResponseSize(2), nSamples);
       
                
                temporalIRCenter = mosaic.temporalImpulseResponseCenterRGB{rgbIndex};
                temporalIRSurround = mosaic.temporalImpulseResponseSurroundRGB{rgbIndex};
                
                fullResponseRSCenter = convn(spResponseCenterRS, temporalIRCenter');
                fullResponseRSSurround = convn(spResponseSurroundRS, temporalIRSurround');
                
                fullResponseRSRGB(:,:,rgbIndex) = fullResponseRSCenter(:,1:nSamples) - fullResponseRSSurround(:,1:nSamples);
                
                %         toc
                
                
                
            end
            
end

% % NEED TO ADD IN EACH RGB CONTRIBUTION
% fullResponse{xcell,ycell,rgbIndex} = reshape(fullResponseRS, spResponseSize(1), spResponseSize(2), size(fullResponseRS,2));
% nlResponse{xcell,ycell,rgbIndex} = exp(mean(fullResponseRS,1));

% fullResponseRS = mean(fullResponseRSRGB,3);
fullResponseRS = sum(fullResponseRSRGB,3);

% fullResponse{xcell,ycell} = reshape(fullResponseRS, spResponseSize(1), spResponseSize(2), size(fullResponseRS,2));
fullResponse{xcell,ycell} = mean(fullResponseRS);
% nlResponse{xcell,ycell} = exp(mean(fullResponseRS,1));

if ~isa(mosaic, 'rgcMosaicLinear')
    
    genFunction = mosaicGet(mosaic, 'generatorFunction');
    nlResponse{xcell,ycell} = genFunction(mean(fullResponseRS,1));
else
    nlResponse{xcell,ycell} = [];
end

end
toc
end