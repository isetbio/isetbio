function  spResponse = spConvolve(mosaic, sptempStimulus)
% spConvolve: a util function of the @rgc parent class, for a separable
% STRF finds the 2D convolution of the spatial RF for every time frame.
% 
% Inputs:
% 
% Outputs:
% 
% Example:
% 
% (c) isetbio
% 09/2015 JRG

stimSize = size(sptempStimulus(:,:,1));
nSamples = size(sptempStimulus,3);

nCells = size(mosaic.spatialRFArray);

% rfSize = size(mosaic.spatialRFArray{1,1});
rfSize = floor(mosaic.receptiveFieldDiameter1STD*ones(2,1));
tic
for xcell = 1:nCells(1)
    for ycell = 1:nCells(2)
        
        for samp = 1:nSamples
            spRF = mosaic.spatialRFArray{xcell,ycell};
            % if stimSize(1) < rfSize(1) && stimSize(2) < rfSize(2)
                spStim = squeeze(sptempStimulus(:,:,samp));
                
                spRFOneDim = mosaic.spatialRFonedim{xcell,ycell};
            % else % need to fix
%                 stimCenterCoords = mosaic.cellCenterLocations{xcell,ycell};
%                 stimX = ceil(stimCenterCoords(1) - floor(rfSize(1)/2)):ceil(stimCenterCoords(1) + floor(rfSize(1)/2));
%                 stimY = ceil(stimCenterCoords(2) - floor(rfSize(2)/2)):ceil(stimCenterCoords(2) + floor(rfSize(2)/2));
%                 spStim = sptempStimulus(stimX,stimY,samp);
            % end
            
            % NEED TO FIX THIS!
            % spResponse{xcell,ycell}(:,:,samp) = conv2(spRF, spStim);
            
            spResponse{xcell,ycell}(:,:,samp) = conv2(spRF, spStim);
            
            
%             tic
%             spOneDim1 = convn(spRFOneDim(1,:), spStim);
%             spOneDim2 = convn(spRFOneDim(2,:), spStim');
%             
% %             spResponse{xcell,ycell}(:,:,samp) = reshape(spOneDim1'*spOneDim2,stimSize);
%             toc
            
        end
    end
end

toc