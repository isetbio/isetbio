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
nSamples = size(spResponse{1,1},3);

nCells = size(mosaic.spatialRFArray);

rfSize = size(mosaic.spatialRFArray{1,1});

% fullResponse = cell(nCells); nlResponse = cell(nCells);
tic

%         tempIRmatrix = convmtx(temporalIR,size(spResponseRS,2));
for xcell = 1:nCells(1)
    for ycell = 1:nCells(2)
        spResponseRS = reshape(spResponse{xcell,ycell}, spResponseSize(1)*spResponseSize(2), nSamples);
        
        temporalIR = mosaic.temporalImpulseResponse;
%         tic
        fullResponseRS = convn(spResponseRS, temporalIR');
%         toc
        
%         tic
%         tempIRmatrix = convmtx(temporalIR,size(spResponseRS,2));
%         fullResponseRS = (tempIRmatrix*spResponseRS');
%         toc
        
        fullResponse{xcell,ycell} = reshape(fullResponseRS, spResponseSize(1), spResponseSize(2), size(fullResponseRS,2));

        nlResponse{xcell,ycell} = exp(mean(fullResponseRS,1));
        
    end
end
toc