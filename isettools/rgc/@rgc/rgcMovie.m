function rgcMovie(obj, outersegment, varargin)

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
%% Movie of RGC response

figure;
set(gcf,'position',[1000  540 893  798]);
cmap = parula(8);



for cellTypeInd = 1:length(obj.mosaic)
    
mosaicall{cellTypeInd} = zeros(245,215,229);

    subplot(3,2,cellTypeInd);
    
    fillIndicesMosaic = rfFill(obj.mosaic{cellTypeInd});
    
    nCells = size(obj.mosaic{cellTypeInd}.cellLocation);
    for xcell = 1:nCells(1)
        for ycell = 1:nCells(2)
            clear x1 y1
            hold on; 
            
            
%             fillIndices = obj.mosaic{cellTypeInd}.spatialRFFill{xcell,ycell};
            fillIndices = fillIndicesMosaic{xcell,ycell};
%             imagesc(obj.mosaic{cellTypeInd}.spatialRFArray{xcell,ycell}(fillIndices));
% 
            [x1 y1] = ind2sub(size(obj.mosaic{cellTypeInd}.sRFcenter{xcell,ycell}),fillIndices);
            x1 = x1 + obj.mosaic{cellTypeInd}.cellLocation{xcell,ycell}(1);
            y1 = y1 + obj.mosaic{cellTypeInd}.cellLocation{xcell,ycell}(2);
%             indall = sub2ind([500 500],x1,y1);
%             mosaicall(round(indall)) = 1;
            indall = ([round(x1(:)) round(y1(:))])';
%             mosaicall(indall(:,1:100)) = 1;
%             mosaicall(fillIndices) = 1;


%                     mosaicall(round(x1(:)),round(y1(1))) = 1;
            
            k = 100:300;%length(obj.mosaic{cellTypeInd}.nlResponse{1,1});

            for ix = 1:length(x1)
%                 for iy = 1%:length(y1)
%                     mosaicall(round(x1(ix)),round(y1(ix)),k) = 1;

                    
%                     mosaicall{cellTypeInd}(round(x1(ix)),round(y1(ix)),k) = log(obj.mosaic{cellTypeInd}.nlResponse{xcell,ycell}(k));

%                     mosaicall{cellTypeInd}(round(x1(ix)),round(y1(ix)),k) = obj.mosaic{cellTypeInd}.nlResponse{xcell,ycell}(k);
% 
%                     mosaicall{cellTypeInd}(round(x1(ix)),round(y1(ix)),ceil(obj.mosaic{cellTypeInd}.spikeResponse{xcell,ycell})) = 1;                    

                    px = ceil(-3*obj.mosaic{cellTypeInd}.rfDiameter + (x1(ix)));
                    py = ceil(-3*obj.mosaic{cellTypeInd}.rfDiameter + (y1(ix)));

                    
                    if isa(obj,'rgcLinear')
                        
                    mosaicall{cellTypeInd}(px,py,k) = log(obj.mosaic{cellTypeInd}.linearResponse{xcell,ycell}(k));
                    else
                        
%                     mosaicall{cellTypeInd}(px,py,ceil(obj.mosaic{cellTypeInd}.spikeResponse{xcell,ycell})) = 1;
                                        
%                     mosaicall{cellTypeInd}(px,py,k) = log(obj.mosaic{cellTypeInd}.nlResponse{xcell,ycell}(k));

                        mosaicall{cellTypeInd}(px,py,k) = (obj.mosaic{cellTypeInd}.psthResponse{xcell,ycell}(k));
                    end
%                 end
            end
            
%             xc = obj.mosaic{cellTypeInd}.cellCenterLocations{xcell,ycell}(1);
%             yc = obj.mosaic{cellTypeInd}.cellCenterLocations{xcell,ycell}(2);
%             mosaicall(fillIndices+sub2ind([500 500],round(xc),round(yc))) = 1;
            
%             [x1 y1] = ind2sub(size(obj.mosaic{cellTypeInd}.spatialRFArray{xcell,ycell}),fillIndices);
%             scatter(x1 + obj.mosaic{cellTypeInd}.cellCenterLocations{xcell,ycell}(1),y1 + obj.mosaic{cellTypeInd}.cellCenterLocations{xcell,ycell}(2));
        end
    end
%     colormap gray
%     figure;
%     imagesc(mosaicall{3}(:,:,140)); drawnow
%     axis equal
%     title(sprintf('%s',obj.mosaic{cellTypeInd}.nameCellType),'fontsize',16);
%     xlabel(sprintf('Distance (\\mum)'),'fontsize',16);
%     ylabel(sprintf('Distance (\\mum)'),'fontsize',16);
end
close;

%%
% h1 = figure;
% set(gcf,'position',[548   606   893   739]);
% vObj = VideoWriter('test.mp4','MPEG-4');
% vObj.FrameRate = 30;
% open(vObj);
% for j = 150:200
%     image(mosaicall{cellTypeInd}(:,:,j)');
% 
%     pr4 = horzcat(obj.mosaic{cellTypeInd}.spatialRFcontours{:,:,1});
%     hold on;
%     plot(pr4(1,:),pr4(2,:),'r','linewidth',2)
% 
%     drawnow
%     F = getframe(h1);
%     writeVideo(vObj,F);
% end
% 
% close(vObj);
% ph=1;
%%
% figure;
tic
h1 = figure;
% set(h1,'position',[0.1 0.1 0.4 0.7])
% set(gcf,'position',[1000  540 893  798]);


set(gcf,'position',[548   606   993   839]);
% vObj = VideoWriter('new2.mj2', 'Archival');

vObj = VideoWriter('testB3.mp4','MPEG-4');
vObj.FrameRate = 30;
vObj.Quality = 100;
 open(vObj);


 sceneRGB = (osGet(outersegment,'rgbData'));
 
%  for cellTypeInd = 1:obj.numberCellTypes
%     pr4 = horzcat(rgc1.mosaic{cellTypeInd}.spatialRFcontours{:});

 plotOrder = [1 4 2 5 3];
% for k = 10


for cellTypeInd = 1:length(obj.mosaic)
    spatialRFcontours{:,:,:,cellTypeInd} = plotContours(obj.mosaic{cellTypeInd});
end

for k = 100:289%length(obj.mosaic{cellTypeInd}.nlResponse{1,1});
    
    fprintf('\b\b\b%02d%%', round(100*k/300));

for cellTypeInd = 1:length(obj.mosaic)
    nCells = size(obj.mosaic{cellTypeInd}.cellLocation);

    
    subplot(2,3,plotOrder(cellTypeInd));
    image(mosaicall{cellTypeInd}(:,:,k)');
    
    % pr4 = horzcat(obj.mosaic{cellTypeInd}.spatialRFcontours{:,:,1});
    
    % spatialRFcontours = plotContours(obj.mosaic{cellTypeInd});
    spatialRFcontoursMosaic = spatialRFcontours{:,:,1,cellTypeInd};
    pr4 = horzcat(spatialRFcontoursMosaic{:,:,1});
    
    hold on;
    plot(pr4(1,:),pr4(2,:),'r','linewidth',2)

    axis equal; axis off;
    title(sprintf('%s',obj.mosaic{cellTypeInd}.cellType),'fontsize',16);

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


