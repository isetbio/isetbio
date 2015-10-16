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

spResponseSize = size(spResponse{1,1}(:,:,1,1));
nSamples = size(spResponse{1,1},3);
channelSize = size(spResponse{1,1},4);

nCells = size(mosaic.cellLocation);

rfSize = size(mosaic.sRFcenter{1,1});

% fullResponse = cell(nCells); nlResponse = cell(nCells);
fprintf('Temporal Convolution, %s:     \n', mosaic.cellType);


%         tempIRmatrix = convmtx(temporalIR,size(spResponseRS,2));


%     fprintf('RGB = %d     \n', rgbIndex);
    tic
for xcell = 1:nCells(1)
    for ycell = 1:nCells(2) 
        
for rgbIndex = 1:channelSize


            if 0%mosaic.temporalImpulseResponseCenterRGB == mosaic.temporalImpulseResponseSurroundRGB
                spResponseRS = reshape(squeeze(spResponse{xcell,ycell}(:,:,:,rgbIndex)), spResponseSize(1)*spResponseSize(2), nSamples);
            
                temporalIR = mosaic.tCenter{rgbIndex};
                
%                  tic
%                 tempIRmatrix = convmtx(temporalIR,size(spResponseRS,2));
%                 fullResponseRS = (tempIRmatrix*spResponseRS');
%                 toc
                
%                 tic
                fullResponseRSRGB(:,:,rgbIndex) = convn(spResponseRS, temporalIR');
%                 toc
                
               
            elseif 1
                
                %         tic
                
                spResponseCenterRS = reshape(squeeze(spResponse{xcell,ycell,1}(:,:,:,rgbIndex)), spResponseSize(1)*spResponseSize(2), nSamples);
                spResponseSurroundRS = reshape(squeeze(spResponse{xcell,ycell,2}(:,:,:,rgbIndex)), spResponseSize(1)*spResponseSize(2), nSamples);
       
                
                temporalIRCenter = mosaic.tCenter{rgbIndex};
                temporalIRSurround = mosaic.tSurround{rgbIndex};
                
                % if strcmpi(mosaic.input, 'rgb') % assume stimulus referred
                
                    fullResponseRSCenter = convn(spResponseCenterRS, temporalIRCenter','same');
                    fullResponseRSSurround = convn(spResponseSurroundRS, temporalIRSurround','same');
                % else % assume cone current referred, do not do any temporal filtering
                
                %     fullResponseRSCenter = spResponseCenterRS;
                %     fullResponseRSSurround = spResponseSurroundRS;
                
                % end
                
                % fullResponseRSRGB(:,:,rgbIndex) = fullResponseRSCenter(:,1:nSamples) - fullResponseRSSurround(:,1:nSamples);
%                 startPoint = length(temporalIRCenter)-1; endPoint = nSamples+length(temporalIRCenter)-1;
%                 fullResponseRSRGB(:,:,rgbIndex) = fullResponseRSCenter(:,startPoint:endPoint) - fullResponseRSSurround(:,startPoint:endPoint);
                fullResponseRSRGB(:,:,rgbIndex) = fullResponseRSCenter - fullResponseRSSurround;
                
                %         toc
                
                
                
            end
            
end

% % NEED TO ADD IN EACH RGB CONTRIBUTION
% fullResponse{xcell,ycell,rgbIndex} = reshape(fullResponseRS, spResponseSize(1), spResponseSize(2), size(fullResponseRS,2));
% nlResponse{xcell,ycell,rgbIndex} = exp(mean(fullResponseRS,1));

% fullResponseRS = mean(fullResponseRSRGB,3);

if ~isa(mosaic,'rgcMosaicSubunit');

    fullResponseRS = sum(fullResponseRSRGB,3);
    
else
    genFunction = mosaicGet(mosaic, 'generatorFunction');
    fullResponseRS = sum(genFunction(fullResponseRSRGB),3);
end

fullResponse{xcell,ycell,1} = mean(fullResponseRS);

% % % fullResponse for RGB
fullResponse{xcell,ycell,2} =  reshape(fullResponseRS, spResponseSize(1), spResponseSize(2), size(fullResponseRS,2));
% % % fullResponse for just G
% fullResponse{xcell,ycell,2} =  reshape(squeeze(fullResponseRSRGB(:,:,2)), spResponseSize(1), spResponseSize(2), size(fullResponseRS,2));
% nlResponse{xcell,ycell} = exp(mean(fullResponseRS,1));

if isa(mosaic, 'rgcMosaicLinear')
    
    nlResponse{xcell,ycell} = [];
    

else
    if isa(mosaic,'rgcMosaicSubunit')
        nlResponse{xcell,ycell} = (mean(fullResponseRS,1));
    else
        genFunction = mosaicGet(mosaic, 'generatorFunction');
        nlResponse{xcell,ycell} = genFunction(mean(fullResponseRS,1));
        
    end
    
end

end
% toc
end
toc