function  [spResponseCenter, spResponseSurround] = spConvolve(mosaic, sptempStimulus, sptempStimulusSurround)
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
% figure;
for rgbIndex = 1:channelSize 
    tic
    
    for xcell = 1:nCells(1)
        for ycell = 1:nCells(2)
            spResponseCenter{xcell,ycell} = zeros([size(mosaic.sRFcenter{xcell,ycell}) nSamples channelSize]);%conv2(spRFcenter, spStim, 'same');
            spResponseSurround{xcell,ycell} = zeros([size(mosaic.sRFcenter{xcell,ycell}) nSamples channelSize]);%conv2(spRFcenter, spStim, 'same');
            for samp = 1:nSamples
                
                % Get RF 2D images
                spRFcenter = mosaic.sRFcenter{xcell,ycell};
                spRFsurround = mosaic.sRFsurround{xcell,ycell};
                
                % Get cell location
                stimCenterCoords = mosaic.cellLocation{xcell,ycell};
                
                % Find the spatial extent of the RF in terms of multiples of rfDiameter
                if isa(mosaic, 'rgcPhys')
                    extent = round(size(mosaic.sRFcenter{1,1},1)/mosaic.rfDiameter);
                    offset = [0 0];%mosaic.cellLocation{1,1};% - floor((extent/2)*mosaic.rfDiameter);
                    
                    %%%% Should rfDimater be size RF center? Yes!
                    stimX =  ceil((stimCenterCoords(1) - floor((extent/2)*mosaic.rfDiameter))/1):floor((stimCenterCoords(1) + floor((extent/2)*mosaic.rfDiameter ))/1);%
                    stimY =  ceil((stimCenterCoords(2) - floor((extent/2)*mosaic.rfDiameter))/1):floor((stimCenterCoords(2) + floor((extent/2)*mosaic.rfDiameter ))/1);%
                
                
                elseif isa(mosaic, 'rgcSubunit')
                    
                    extent = 1;%round(size(mosaic.sRFcenter{1,1},1)/mosaic.rfDiameter);
                    %                     offset = mosaic.cellLocation{1,1};% - floor((extent/2)*mosaic.rfDiameter);
                    if mosaic.cellLocation{1,1}(1) >0
                        offset(1) = ceil(mosaic.cellLocation{1,1}(1));
                    else
                        offset(1) = floor(mosaic.cellLocation{1,1}(1));
                    end
                    
                    if mosaic.cellLocation{1,1}(2) >0
                        offset(2) = ceil(mosaic.cellLocation{1,1}(2));
                    else
                        offset(2) = floor(mosaic.cellLocation{1,1}(2));
                    end
%                     offset = offset+2*(rand(1,2)-1);
                    %%%% Should rfDimater be size RF center? Yes!
                    stimX =  floor((stimCenterCoords(1) - floor((extent/2)*size(mosaic.sRFcenter{1,1},1)))/1):floor((stimCenterCoords(1) + floor((extent/2)*size(mosaic.sRFcenter{1,1},1) ))/1);%
                    stimY =  floor((stimCenterCoords(2) - floor((extent/2)*size(mosaic.sRFcenter{1,1},2)))/1):floor((stimCenterCoords(2) + floor((extent/2)*size(mosaic.sRFcenter{1,1},2)))/1);%
                    
                else
                    extent = 1;%2.5/2;%round(size(mosaic.sRFcenter{1,1},1)/mosaic.rfDiameter);
                    offset = mosaic.cellLocation{1,1};% - floor((extent/2)*mosaic.rfDiameter);
                    
                    %%%% Should rfDimater be size RF center? Yes!
%                     stimX =  ceil((stimCenterCoords(1) - floor((extent/2)*mosaic.rfDiameter))/1):floor((stimCenterCoords(1) + floor((extent/2)*mosaic.rfDiameter ))/1);%
%                     stimY =  ceil((stimCenterCoords(2) - floor((extent/2)*mosaic.rfDiameter))/1):floor((stimCenterCoords(2) + floor((extent/2)*mosaic.rfDiameter ))/1);%
                    stimX =  ceil((stimCenterCoords(1) - floor((extent/2)*size(mosaic.sRFcenter{1,1},1)))/1):floor((stimCenterCoords(1) + floor((extent/2)*size(mosaic.sRFcenter{1,1},1) ))/1);%
                    stimY =  ceil((stimCenterCoords(2) - floor((extent/2)*size(mosaic.sRFcenter{1,1},2)))/1):floor((stimCenterCoords(2) + floor((extent/2)*size(mosaic.sRFcenter{1,1},2)))/1);%
                
                end
                
                % Extract x and y coordinates within spatial extent offset by the center coordinate
                % make non-rectangular? follow exact RF contours?
                

                % Ensure indices are within size of stimulus
%                 gz = find(stimX>=1 & stimY>=1 & stimX<=size(sptempStimulus,1) & stimY<=size(sptempStimulus,2) );
                
                gz = find((stimX-offset(1))>=1 & (stimY-offset(2))>=1 & (stimX-offset(1))<=size(sptempStimulus,1) & (stimY-offset(2))<=size(sptempStimulus,2) );
                gzout(xcell,ycell,1:length(gz(:))) = gz(:);
                % Extract 2D image
%                 spStim = squeeze(sptempStimulus(stimX(gz),stimY(gz),samp,rgbIndex));
                spStim = squeeze(sptempStimulus(floor(stimX(gz)-offset(1)),floor(stimY(gz)-offset(2)),samp,rgbIndex));
                % spStim = spStim/max(spStim(:)); spStim = spStim - mean(spStim(:));
                            
                if isa(mosaic, 'rgcSubunit') && exist('spStimSurr','var')
                    
                    spStimSurr = squeeze(sptempStimulusSurround(floor(stimX(gz)-offset(1)),floor(stimY(gz)-offset(2)),samp,rgbIndex));
                end
                
%                 [m1,m2] = meshgrid(floor(stimX(gz)-offset(1)),floor(stimY(gz)-offset(2)));
%                 hold on; scatter(m1(:)+2*(rand-.5),m2(:)+2*(rand-.5),'x');
                % Convolve for a single temporal frame
                
                if isa(mosaic, 'rgcPhys')
%                     spResponseCenter{xcell,ycell}(:,:,samp,rgbIndex) = zeros(size(spStim));%conv2(spRFcenter, spStim, 'same');
                    spResponseCenter{xcell,ycell}(gz,gz,samp,rgbIndex) = spRFcenter(gz,gz).*spStim;%conv2(spRFcenter, spStim, 'same');
%                     spResponseSurround{xcell,ycell}(gz,gz,samp,rgbIndex) = zeros(size(spStim));%conv2(spRFsurround, spStim, 'same');
                
                elseif isa(mosaic, 'rgcSubunit')
%                     
%                     spRC = conv2(spRFcenter, spStim-1*mean(spStim(:)), 'same');
%                     spRS = conv2(spRFsurround, spStim-1*mean(spStim(:)), 'same');

%                     spResponseCenter{xcell,ycell}(gz,gz,samp,rgbIndex) = mosaic.rectifyFunction(sum(spRC(:)));
%                     spResponseSurround{xcell,ycell}(gz,gz,samp,rgbIndex) = mosaic.rectifyFunction(sum(spRS(:)));

                    spRC = (spRFcenter(gz,gz).*(spStim-0*mean(spStim(:))));
                    
                    if exist('spStimSurr','var')
                        spRS = (spRFsurround(gz,gz).*-0*(spStimSurr-1*mean(spStimSurr(:))));
                        spResponseCenter{xcell,ycell}(gz,gz,samp,rgbIndex) = mosaic.rectifyFunction(sum(spRC(:)))./length(gz)^2;
                        spResponseSurround{xcell,ycell}(gz,gz,samp,rgbIndex) = mosaic.rectifyFunction(sum(spRS(:)))./length(gz)^2;

                    else 
                        spRS = (spRFsurround(gz,gz).*(spStim-0*mean(spStim(:))));
                        
%                         spResponseCenter{xcell,ycell}(gz,gz,samp,rgbIndex) = abs(sum(spRC(:)))./1;%length(gz)^2;
%                         spResponseSurround{xcell,ycell}(gz,gz,samp,rgbIndex) = abs(sum(spRS(:)))./1;%length(gz)^2;

%                       spResponseCenter{xcell,ycell}(1,1,samp,rgbIndex) = 10*sum(mosaic.rectifyFunction(spRC(:)))./length(gz)^2;
%                         spResponseSurround{xcell,ycell}(1,1,samp,rgbIndex) = 10.*sum(mosaic.rectifyFunction(spRS(:)))./length(gz)^2;

                        spResponseCenter{xcell,ycell}(gz,gz,samp,rgbIndex) = 10*(mosaic.rectifyFunction(spRC))./length(gz)^2;
                        spResponseSurround{xcell,ycell}(gz,gz,samp,rgbIndex) = 10*(mosaic.rectifyFunction(spRS))./length(gz)^2;

                    end
                    

                else
                    
%                     if xcell == 5

%                     convCenter = conv2(spRFcenter, spStim, 'full');
%                     convSurround = conv2(spRFsurround, spStim, 'full');
%                     sizeCC = size(convCenter); sizeStim = size(spStim);
                    
%                     spResponseCenter{xcell,ycell}(:,:,samp,rgbIndex) = zeros(sizeStim);
%                     spResponseSurround{xcell,ycell}(:,:,samp,rgbIndex) = zeros(sizeStim);
                    
                    
                    
%                     spResponseCenter{xcell,ycell}(ceil(sizeCC(1)/2)-ceil(sizeStim(1)/2):floor(sizeCC(1)/2)+floor(sizeStim(1)/2),ceil(sizeCC(2)/2)-ceil(sizeStim(2)/2):floor(sizeCC(2)/2)+floor(sizeStim(2)/2),samp,rgbIndex) = convCenter(...
%                         ceil(sizeCC(1)/2)-ceil(sizeStim(1)/2):floor(sizeCC(1)/2)+floor(sizeStim(1)/2),...
%                         ceil(sizeCC(2)/2)-ceil(sizeStim(2)/2):floor(sizeCC(2)/2)+floor(sizeStim(2)/2));
%                     
%                     spResponseSurround{xcell,ycell}(ceil(sizeCC(1)/2)-ceil(sizeStim(1)/2):floor(sizeCC(1)/2)+floor(sizeStim(1)/2),ceil(sizeCC(2)/2)-ceil(sizeStim(2)/2):floor(sizeCC(2)/2)+floor(sizeStim(2)/2),samp,rgbIndex) = convSurround(...
%                         ceil(sizeCC(1)/2)-ceil(sizeStim(1)/2):floor(sizeCC(1)/2)+floor(sizeStim(1)/2),...
%                         ceil(sizeCC(2)/2)-ceil(sizeStim(2)/2):floor(sizeCC(2)/2)+floor(sizeStim(2)/2));
                    
                    
%                     sizeResp = size(spResponseCenter{xcell,ycell}(:,:,samp,rgbIndex) );
%                     
%                     if (sizeResp(1) ~= sizeStim(1)) && (sizeResp(2) ~= sizeStim(2))
                        
                    
%                     spResponseSurround{xcell,ycell}(:,:,samp,rgbIndex) = conv2(spRFsurround, spStim, 'same');


% % % % % % % % % % % % % % % % % 

%                     spResponseCenter{xcell,ycell}(:,:,samp,rgbIndex) = conv2(spRFcenter, spStim-1*mean(spStim(:)), 'same');
%                     spResponseSurround{xcell,ycell}(:,:,samp,rgbIndex) = conv2(spRFsurround, spStim-1*mean(spStim(:)), 'same');
% 
                    spResponseCenter{xcell,ycell}(gz,gz,samp,rgbIndex) = (spRFcenter(gz,gz).*(spStim-0*mean(spStim(:))));
                    spResponseSurround{xcell,ycell}(gz,gz,samp,rgbIndex) = (spRFsurround(gz,gz).*(spStim-0*mean(spStim(:))));
 
% % % % % % % % % % % % %                     
                    
                    
                    %                     else
%                              spResponseCenter{xcell,ycell}(:,:,samp,rgbIndex) = zeros(size(conv2(spRFcenter, spStim, 'same')));
%                     spResponseSurround{xcell,ycell}(:,:,samp,rgbIndex) =zeros(size(conv2(spRFsurround, spStim, 'same')));
%                     end
                    
%                     spResponseCenter{xcell,ycell}(:,:,samp,rgbIndex) = spRFcenter(gz,gz).*spStim;
%                     spResponseSurround{xcell,ycell}(:,:,samp,rgbIndex) = spRFsurround(gz,gz).*spStim;

                end
                clear  gz
            end
        end
    end
    
    toc
end
% figure; hold on
% for i1 = 1:5
%     for j1=1:5
%         [m1,m2] = meshgrid(gzout(i1,j1,:) + mosaic.cellLocation{i1,j1}(1),gzout(i1,j1,:)+ mosaic.cellLocation{i1,j1}(2));
%         scatter(m1(:),m2(:),'x');
%     end
% end