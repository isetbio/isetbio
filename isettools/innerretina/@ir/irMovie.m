function irMovie(obj, outersegment, varargin)
% 
% rgcMovie: a method of @rgc that generates a movie of the stimulus and the
% responses of the rgc mosaics. The type of response shown in the movie is
% determined by user input; the default is the PSTH, which represents
% the output of the neurons with a continuous value over time. 
% 
% Inputs: the rgc object, the outersegment object (to show the stimulus)
% TO DO: make flexible for osLinear input
% 
% Outputs: a *.mp4 file containing the movie of the output responses.
% TO DO: allow the user to input a name for the movie file.
% 
% Example:
%   rgcMovie(rgc1, os);
% 
% (c) isetbio
% 09/2015 JRG

%% Build movies for each mosaic
% Since the rgc object contains at least 5 mosaics, the movie for each
% mosaic must be made before they can be integrated into the same figure.
% Alternatively, this could be done with something like ffmpeg.

figure;
set(gcf,'position',[1000  540 893  798]);
% cmap = parula(8);

subsamp = 1;

for cellTypeInd = 1:length(obj.mosaic)
        
    subplot(3,2,cellTypeInd);
    
    mosaicall{cellTypeInd} = zeros(100,100,300);
    
    % rfFill is a utility function that colors in cells in the mosaic
    fillIndicesMosaic = rfFill(obj.mosaic{cellTypeInd});
    
    nCells = size(obj.mosaic{cellTypeInd}.cellLocation);
        
    extent = .5*round(size(obj.mosaic{cellTypeInd}.sRFcenter{1,1},1)/obj.mosaic{cellTypeInd}.rfDiameter);
                
    for xcell = 1:nCells(1)
        for ycell = 1:nCells(2)
            clear x1 y1
            hold on; 
            
            fillIndices = fillIndicesMosaic{xcell,ycell};
 
            % Shift fill indices to cell center coordinate
            [x1 y1] = ind2sub(size(obj.mosaic{cellTypeInd}.sRFcenter{xcell,ycell}),fillIndices);
            x1 = x1 + obj.mosaic{cellTypeInd}.cellLocation{xcell,ycell}(1);
            y1 = y1 + obj.mosaic{cellTypeInd}.cellLocation{xcell,ycell}(2);
            
            k = 1:subsamp:500;%length(obj.mosaic{cellTypeInd}.nlResponse{1,1});

            % Color fill indices by magntidue of response
            % Need to loop instead of vectorize for some unknown reason
            for ix = 1:length(x1)

                px = max([1,ceil(-extent*obj.mosaic{cellTypeInd}.rfDiameter + (x1(ix)))]);
                py = max([1,ceil(-extent*obj.mosaic{cellTypeInd}.rfDiameter + (y1(ix)))]);
                               
                if isa(obj,'rgcLinear')                 
                    mosaicall{cellTypeInd}(px,py,(k-1)/subsamp+1) = log(obj.mosaic{cellTypeInd}.linearResponse{xcell,ycell}(k));
                else
                    
                    mosaicall{cellTypeInd}(px,py,ceil(10*obj.mosaic{cellTypeInd}.spikeResponse{xcell,ycell})) = 100;                   
                    % mosaicall{cellTypeInd}(px,py,k) = log(obj.mosaic{cellTypeInd}.nlResponse{xcell,ycell}(k));              
                    
%                     mosaicall{cellTypeInd}(px,py,(k-1)/subsamp+1) = (obj.mosaic{cellTypeInd}.psthResponse{xcell,ycell}(k));
                         
%                     mosaicall{cellTypeInd}(px,py,ceil(obj.mosaic{cellTypeInd}.spikeResponse{xcell,ycell,1,1})) = 1;
                    % Bring max values to front
                    % m1 = mosaicall{cellTypeInd}(px,py,k); m2 = (obj.mosaic{cellTypeInd}.psthResponse{xcell,ycell}(k));
                    % mosaicall{cellTypeInd}(px,py,k) = max([m1(:)'; m2(:)'  ]);
                    % clear m1 m2
                end
            end
            
        end
    end
end
close;


%% Build figure with RF contours and subplot and capture using gcf

% figure;
tic
h1 = figure;
% set(h1,'position',[0.1 0.1 0.4 0.7])
% set(gcf,'position',[1000  540 893  798]);
set(gcf,'position',[548   606   993   839]);

% Initialize video file
vObj = VideoWriter('barJan28.mp4','MPEG-4');
vObj.FrameRate = 30;
vObj.Quality = 100;
open(vObj);

if isa(outersegment, 'osIdentity')
    sceneRGB = (osGet(outersegment,'rgbData'));
elseif strcmpi(outersegment.type, 'sensor')
    sceneRGB = sensorGet(outersegment,'volts');
end
 
% Rearrange mosaic types for movie display
plotOrder = [1 4 2 5 3];

% Generate spatial RF contours
for cellTypeInd = 1:length(obj.mosaic)
    spatialRFcontours{:,:,:,cellTypeInd} = plotContours(obj.mosaic{cellTypeInd});
end

% Set axes limits
xAxisLimit = round(1.33*size(squeeze(sceneRGB(:,:,1,1)),1));
yAxisLimit = round(1.33*size(squeeze(sceneRGB(:,:,1,1)),2));

% Build each frame and gcf
for k = 1:subsamp:500%length(obj.mosaic{cellTypeInd}.nlResponse{1,1});
    
    fprintf('\b\b\b%02d%%', round(100*k/300));
    
    for cellTypeInd = 1:length(obj.mosaic)
        nCells = size(obj.mosaic{cellTypeInd}.cellLocation);
               
        subplot(2,3,plotOrder(cellTypeInd));
        
        % Draw fill using image
        image(mosaicall{cellTypeInd}(:,:,(k-1)/subsamp+1)'); colormap gray
        
        % Draw RF contours on same image
%         plot(spatialRFcontours{1,1,1}(1,2:end),spatialRFcontours{xcell,ycell,1}(2,2:end))
%             spatialRFcontours{xcell,ycell,1}(2,2:end),...
%             'color',cmap(cellTypeInd,:));
        spatialRFcontoursMosaic = spatialRFcontours{:,:,1,cellTypeInd};
        spatialRFcontoursMosaicArr = horzcat(spatialRFcontoursMosaic{:,:,1});
        
        hold on;
        plot(spatialRFcontoursMosaicArr(1,:),spatialRFcontoursMosaicArr(2,:),'r','linewidth',2)
        
        axis equal; axis off; axis([0 xAxisLimit 0 yAxisLimit]);
        title(sprintf('%s',obj.mosaic{cellTypeInd}.cellType),'fontsize',16);
        
        % Draw frame from stimulus movie
%         subplot(2,3,6); image(squeeze(sceneRGB(:,:,1+floor((k-0)/10),:))); axis equal; axis off;
        subplot(2,3,6); image(squeeze(sceneRGB(:,:,1+floor((k-0)/10),:))); axis equal; axis off;
        title('Stimulus', 'fontsize', 16);
             
    end

    drawnow
    F = getframe(h1);
    writeVideo(vObj,F);
    
end
close(vObj)
toc;
fprintf('     \n');
ph = 1;


