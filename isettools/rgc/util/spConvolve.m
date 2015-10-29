function  [spResponseCenter, spResponseSurround] = spConvolve(mosaic, sptempStimulus)
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

nSamples = size(sptempStimulus,3);
channelSize = size(sptempStimulus,4);

nCells = size(mosaic.cellLocation);

fprintf('     \n');
fprintf('Spatial Convolution, %s:     \n', mosaic.cellType);

for rgbIndex = 1:channelSize
    tic
    
    for xcell = 1:nCells(1)
        for ycell = 1:nCells(2)
            
            for samp = 1:nSamples
                
                spRFcenter = mosaic.sRFcenter{xcell,ycell};
                spRFsurround = mosaic.sRFsurround{xcell,ycell};
                
                stimCenterCoords = mosaic.cellLocation{xcell,ycell};
                extent = round(size(mosaic.sRFcenter{1,1},1)/mosaic.rfDiameter);
                
                % make non-rectangular? follow exact RF contours?
                stimX =  ceil((stimCenterCoords(1) - floor((extent/2)*mosaic.rfDiameter))/1):floor((stimCenterCoords(1) + floor((extent/2)*mosaic.rfDiameter ))/1);%
                stimY =  ceil((stimCenterCoords(2) - floor((extent/2)*mosaic.rfDiameter))/1):floor((stimCenterCoords(2) + floor((extent/2)*mosaic.rfDiameter ))/1);%
                                
                gz = find(stimX>=1 & stimY>=1 & stimX<=size(sptempStimulus,1) & stimY<=size(sptempStimulus,2) );
                spStim = squeeze(sptempStimulus(stimX(gz),stimY(gz),samp,rgbIndex));
                                                       
                spResponseCenter{xcell,ycell}(:,:,samp,rgbIndex) = conv2(spRFcenter, spStim, 'same');
                spResponseSurround{xcell,ycell}(:,:,samp,rgbIndex) = conv2(spRFsurround, spStim, 'same');
                
                clear  gz
            end
        end
    end
    
    toc
end