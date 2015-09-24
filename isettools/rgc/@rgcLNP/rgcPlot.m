function rgcPlot(obj, sensor, outersegment, varargin)
% rgcPlot: a method of @rgcLNP that plots the retinal location of the
% patch, the 1 standard deviation contours of the RGC RF mosaic and a movie
% of the RGC response.
% 
% Inputs:
% 
% Outputs:
% 
% Example:
% 
% (c) isetbio
% 09/2015 JRG


%% Plot retinal location
   
% Plot the TEE and the location of the retinal patch.
plotPatchEccentricity(obj.patchLocationPolarAngleDegrees, obj.patchLocationPolarRadiusMicrometers, obj.eyeLeftOrRight, obj.temporalEquivEcc);

%% Plot the mosaic of each RGC type


figure;
set(gcf,'position',[1000  540 893  798]);
cmap = parula(8);

for cellTypeInd = 1:obj.numberCellTypes
    
    subplot(3,2,cellTypeInd);
    
    nCells = size(obj.mosaic{cellTypeInd}.spatialRFArray);
    for xcell = 1:nCells(1)
        for ycell = 1:nCells(2)
            hold on; 
            plot(obj.mosaic{cellTypeInd}.spatialRFcontours{xcell,ycell}(1,2:end),...
                obj.mosaic{cellTypeInd}.spatialRFcontours{xcell,ycell}(2,2:end),...
                'color',cmap(cellTypeInd,:));
        end
    end
    axis equal
    title(sprintf('%s',obj.mosaic{cellTypeInd}.nameCellType),'fontsize',16);
    xlabel(sprintf('Distance (\\mum)'),'fontsize',16);
    ylabel(sprintf('Distance (\\mum)'),'fontsize',16);
end

%% Movie of RGC response

figure;
set(gcf,'position',[1000  540 893  798]);
cmap = parula(8);

for cellTypeInd = 1:obj.numberCellTypes
    
mosaicall{cellTypeInd} = zeros(245,215,229);

    subplot(3,2,cellTypeInd);
    
    nCells = size(obj.mosaic{cellTypeInd}.spatialRFArray);
    for xcell = 1:nCells(1)
        for ycell = 1:nCells(2)
            clear x1 y1
            hold on; 
            fillIndices = obj.mosaic{cellTypeInd}.spatialRFFill{xcell,ycell};
%             imagesc(obj.mosaic{cellTypeInd}.spatialRFArray{xcell,ycell}(fillIndices));
% 
            [x1 y1] = ind2sub(size(obj.mosaic{cellTypeInd}.spatialRFArray{xcell,ycell}),fillIndices);
            x1 = x1 + obj.mosaic{cellTypeInd}.cellCenterLocations{xcell,ycell}(1);
            y1 = y1 + obj.mosaic{cellTypeInd}.cellCenterLocations{xcell,ycell}(2);
%             indall = sub2ind([500 500],x1,y1);
%             mosaicall(round(indall)) = 1;
            indall = ([round(x1(:)) round(y1(:))])';
%             mosaicall(indall(:,1:100)) = 1;
%             mosaicall(fillIndices) = 1;


%                     mosaicall(round(x1(:)),round(y1(1))) = 1;
            
            k = 1:10;%length(obj.mosaic{cellTypeInd}.nlResponse{1,1});

            for ix = 1:length(x1)
%                 for iy = 1%:length(y1)
%                     mosaicall(round(x1(ix)),round(y1(ix)),k) = 1;

                    
%                     mosaicall{cellTypeInd}(round(x1(ix)),round(y1(ix)),k) = log(obj.mosaic{cellTypeInd}.nlResponse{xcell,ycell}(k));

%                     mosaicall{cellTypeInd}(round(x1(ix)),round(y1(ix)),k) = obj.mosaic{cellTypeInd}.nlResponse{xcell,ycell}(k);
% 
%                     mosaicall{cellTypeInd}(round(x1(ix)),round(y1(ix)),ceil(obj.mosaic{cellTypeInd}.spikeResponse{xcell,ycell})) = 1;                    

                    px = ceil(-3*obj.mosaic{cellTypeInd}.receptiveFieldDiameter1STD + (x1(ix)));
                    py = ceil(-3*obj.mosaic{cellTypeInd}.receptiveFieldDiameter1STD + (y1(ix)));

                    mosaicall{cellTypeInd}(px,py,ceil(obj.mosaic{cellTypeInd}.spikeResponse{xcell,ycell})) = 1;
                                        
%                     mosaicall{cellTypeInd}(px,py,k) = log(obj.mosaic{cellTypeInd}.nlResponse{xcell,ycell}(k));
%                 end
            end
            
%             xc = obj.mosaic{cellTypeInd}.cellCenterLocations{xcell,ycell}(1);
%             yc = obj.mosaic{cellTypeInd}.cellCenterLocations{xcell,ycell}(2);
%             mosaicall(fillIndices+sub2ind([500 500],round(xc),round(yc))) = 1;
            
%             [x1 y1] = ind2sub(size(obj.mosaic{cellTypeInd}.spatialRFArray{xcell,ycell}),fillIndices);
%             scatter(x1 + obj.mosaic{cellTypeInd}.cellCenterLocations{xcell,ycell}(1),y1 + obj.mosaic{cellTypeInd}.cellCenterLocations{xcell,ycell}(2));
        end
    end
    colormap gray
%     figure;
%     imagesc(mosaicall{3}(:,:,140)); drawnow
    axis equal
    title(sprintf('%s',obj.mosaic{cellTypeInd}.nameCellType),'fontsize',16);
    xlabel(sprintf('Distance (\\mum)'),'fontsize',16);
    ylabel(sprintf('Distance (\\mum)'),'fontsize',16);
end

% colormap(gray);
% nframes = size(absorptions,3);
% % Record the movie
% for j = 1:step:nframes 
%     image(absorptions(:,:,j)); drawnow;
%     title('Cone absorptions')
% %     F = getframe;
% %     writegit puVideo(vObj,F);
% end
%%
% figure;
tic
h1 = figure;
% set(h1,'position',[0.1 0.1 0.4 0.7])
% set(gcf,'position',[1000  540 893  798]);


set(gcf,'position',[548   606   893   739]);
% vObj = VideoWriter('new2.mj2', 'Archival');

vObj = VideoWriter('test.mp4','MPEG-4');
vObj.FrameRate = 30;
 open(vObj);


 sceneRGB = (osGet(outersegment,'rgbData'));
 
%  for cellTypeInd = 1:obj.numberCellTypes
%     pr4 = horzcat(rgc1.mosaic{cellTypeInd}.spatialRFcontours{:});

 plotOrder = [1 4 2 5 3];
% for k = 10

for k = 1:10%length(obj.mosaic{cellTypeInd}.nlResponse{1,1});
for cellTypeInd = 1:obj.numberCellTypes
        nCells = size(obj.mosaic{cellTypeInd}.spatialRFArray);

    
    subplot(2,3,plotOrder(cellTypeInd));
    imagesc(mosaicall{cellTypeInd}(:,:,k)');
    
    pr4 = horzcat(obj.mosaic{cellTypeInd}.spatialRFcontours{:});
    hold on;
    plot(pr4(1,:),pr4(2,:),'r','linewidth',2)
    
%         for xcell = 1:nCells(1)
%         for ycell = 1:nCells(2)
%                 hold on; 
% %             plot(3*obj.mosaic{cellTypeInd}.receptiveFieldDiameter1STD+obj.mosaic{cellTypeInd}.spatialRFcontours{xcell,ycell}(1,2:end),...
% %                 3*obj.mosaic{cellTypeInd}.receptiveFieldDiameter1STD+obj.mosaic{cellTypeInd}.spatialRFcontours{xcell,ycell}(2,2:end),...
% %                 'color','r','linewidth',2);
%                 
%             plot(0*obj.mosaic{cellTypeInd}.receptiveFieldDiameter1STD+obj.mosaic{cellTypeInd}.spatialRFcontours{xcell,ycell}(1,2:end),...
%                 0*obj.mosaic{cellTypeInd}.receptiveFieldDiameter1STD+obj.mosaic{cellTypeInd}.spatialRFcontours{xcell,ycell}(2,2:end),...
%                  'color','r','linewidth',2);
%         end
%         end
    
    caxis([0 1]);
    axis equal; axis off;
    title(sprintf('%s',obj.mosaic{cellTypeInd}.nameCellType),'fontsize',16);

    subplot(2,3,6); image(squeeze(sceneRGB(:,:,:,1+mod(k,10)))); axis equal; axis off;
    title('Stimulus', 'fontsize', 16);
    
    drawnow
    F = getframe(h1);
    writeVideo(vObj,F);
end
end
toc
close(vObj)
ph = 1;


