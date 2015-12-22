function  [spResponseCenter, spResponseSurround] = spConvolve(mosaic, sptempStimulus)
% spConvolve: a util function of the @rgc parent class, for a separable
% STRF finds the 2D convolution of the spatial RF for every time frame.
% 
%     [spResponseCenter, spResponseSurround] = spConvolve(mosaic, spTempStim);
%   
% This function extracts the relevant x and y coordinates of each temporal
% frame of the stimulus and convolves that 2D image with the spatial center
% and surround RFs of each cell for all temporal frames.
% 
% Inputs: mosaic object, spTempStim (x, y, t)
% 
% Outputs: spatial response of the 2D center and surround receptive fields.
% 
% Example:
%     [spResponseCenter, spResponseSurround] = spConvolve(rgc1.mosaic{1}, spTempStim);
%     
% (c) isetbio
% 09/2015 JRG

nSamples = size(sptempStimulus,3);
channelSize = size(sptempStimulus,4); % if RGB, channelSize = 3

nCells = size(mosaic.cellLocation);

fprintf('     \n');
fprintf('Spatial Convolution, %s:     \n', mosaic.cellType);

for rgbIndex = 1:channelSize 
    tic
    
    for xcell = 1:nCells(1)
        for ycell = 1:nCells(2)
            
            for samp = 1:nSamples
                
                % Get RF 2D images
                spRFcenter = mosaic.sRFcenter{xcell,ycell};
                spRFsurround = mosaic.sRFsurround{xcell,ycell};
                
                % Get cell location
                stimCenterCoords = mosaic.cellLocation{xcell,ycell};
                
                % Find the spatial extent of the RF in terms of multiples of rfDiameter
                extent = round(size(mosaic.sRFcenter{1,1},1)/mosaic.rfDiameter);
                
                % Extract x and y coordinates within spatial extent offset by the center coordinate
                % make non-rectangular? follow exact RF contours?
                stimX =  ceil((stimCenterCoords(1) - floor((extent/2)*mosaic.rfDiameter))/1):floor((stimCenterCoords(1) + floor((extent/2)*mosaic.rfDiameter ))/1);%
                stimY =  ceil((stimCenterCoords(2) - floor((extent/2)*mosaic.rfDiameter))/1):floor((stimCenterCoords(2) + floor((extent/2)*mosaic.rfDiameter ))/1);%
                
                % Ensure indices are within size of stimulus
                gz = find(stimX>=1 & stimY>=1 & stimX<=size(sptempStimulus,1) & stimY<=size(sptempStimulus,2) );
                
                % Extract 2D image
                spStim = squeeze(sptempStimulus(stimX(gz),stimY(gz),samp,rgbIndex));
                % spStim = spStim/max(spStim(:)); spStim = spStim - mean(spStim(:));
                                          
                % Convolve for a single temporal frame
                spResponseCenter{xcell,ycell}(:,:,samp,rgbIndex) = spRFcenter.*spStim;%conv2(spRFcenter, spStim, 'same');
                spResponseSurround{xcell,ycell}(:,:,samp,rgbIndex) = zeros(size(spStim));%conv2(spRFsurround, spStim, 'same');
                
                clear  gz
            end
        end
    end
    
    toc
end