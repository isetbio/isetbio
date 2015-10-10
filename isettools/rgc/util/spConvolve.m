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

nSamples = size(sptempStimulus,3);
channelSize = size(sptempStimulus,4);

stimSize = size(sptempStimulus(:,:,1));

nCells = size(mosaic.cellLocation);

% rfSize = size(mosaic.spatialRFArray{1,1});
rfSize = floor(mosaic.rfDiameter*ones(2,1));
% tic
fprintf('     \n');
fprintf('Spatial Convolution, %s:     \n', mosaic.cellType);

% sptempStimulus = sptempStimulus - mean(sptempStimulus(:));

for rgbIndex = 1:channelSize
    tic    
%     fprintf('RGB = %d     \n', rgbIndex);
for xcell = 1:nCells(1)
    for ycell = 1:nCells(2)
        
        for samp = 1:nSamples
            % spRF = mosaic.spatialRFArray{xcell,ycell};
            spRFcenter = mosaic.sRFcenter{xcell,ycell};
            spRFsurround = mosaic.sRFsurround{xcell,ycell};
            % if stimSize(1) < rfSize(1) && stimSize(2) < rfSize(2)

%                 spStim = squeeze(sptempStimulus(:,:,samp,rgbIndex));

%                 spStim = spStim - mean(spStim(:));
                
                % spRFOneDim = mosaic.spatialRFonedim{xcell,ycell};
            % else % need to fix
                stimCenterCoords = mosaic.cellLocation{xcell,ycell};
                extent = round(size(mosaic.sRFcenter{1,1},1)/mosaic.rfDiameter);
                stimX =  (ceil(stimCenterCoords(1) - [1; 1]*floor(mosaic.rfDiameter/2)):floor(stimCenterCoords(1) + [1; 1]*floor(mosaic.rfDiameter/2))) - floor((extent/2)*mosaic.rfDiameter);
                
                stimY =  (ceil(stimCenterCoords(2) - [1; 1]*floor(mosaic.rfDiameter/2)):floor(stimCenterCoords(2) + [1; 1]*floor(mosaic.rfDiameter/2))) - floor((extent/2)*mosaic.rfDiameter);
                
                spStim = zeros(length(stimX),length(stimY));
                
                [stimXgrid,stimYgrid] = meshgrid(stimX,stimY);
                
                lz = find(stimXgrid<1|stimYgrid<1|stimXgrid>length(stimX)|stimYgrid>length(stimY));
                if length(lz) > 0
                    spStim(lz) = 0.5;
                    gz = find(stimXgrid>=1 & stimYgrid>=1);
                    % [gzx gzy] = ind2sub([28,28],gz);
                    spTempFull = sptempStimulus(:,:,samp,rgbIndex);
                    spStim(gz) = spTempFull(gz); % squeeze(sptempStimulus(-14+gzx,-14+gzy,samp,rgbIndex));
                else
                    spStim = squeeze(sptempStimulus(stimX,stimY,samp,rgbIndex));
                end
                
%                 stimX = ceil(stimCenterCoords(1) - floor(rfSize(1)/2)):ceil(stimCenterCoords(1) + floor(rfSize(1)/2));
%                 stimY = ceil(stimCenterCoords(2) - floor(rfSize(2)/2)):ceil(stimCenterCoords(2) + floor(rfSize(2)/2));
%                 spStim = sptempStimulus(stimX,stimY,samp);
            % end
            
            % NEED TO FIX THIS!
            % spResponse{xcell,ycell}(:,:,samp) = conv2(spRF, spStim);
            
            if 0
            
                spResponse{xcell,ycell}(:,:,rgbIndex,samp) = conv2(spRF, spStim);
            elseif 1 
                spResponse{xcell,ycell,1}(:,:,samp,rgbIndex) = conv2(spRFcenter, spStim, 'same');            
                spResponse{xcell,ycell,2}(:,:,samp,rgbIndex) = conv2(spRFsurround, spStim, 'same');
                
%                 spResponse{xcell,ycell,1}(:,:,samp,rgbIndex) = conv2(spStim, spRFcenter, 'same');            
%                 spResponse{xcell,ycell,2}(:,:,samp,rgbIndex) = conv2(spStim, spRFsurround, 'same');
            end
%             tic
%             spOneDim1 = convn(spRFOneDim(1,:), spStim);
%             spOneDim2 = convn(spRFOneDim(2,:), spStim');
%             
% %             spResponse{xcell,ycell}(:,:,samp) = reshape(spOneDim1'*spOneDim2,stimSize);

        end
    end
end

toc
end