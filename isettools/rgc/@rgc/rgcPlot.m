function rgcPlot(obj, varargin)
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
   
% % Plot the TEE and the location of the retinal patch.
% plotPatchEccentricity(obj.patchLocationPolarAngleDegrees, obj.patchLocationPolarRadiusMicrometers, obj.eyeLeftOrRight, obj.temporalEquivEcc);

%% Plot the mosaic of each RGC type


figure;
set(gcf,'position',[1000  540 893  798]);
cmap = parula(16);

for cellTypeInd = 1:obj.numberCellTypes
    
    subplot(3,2,cellTypeInd);
    
    nCells = size(obj.mosaic{cellTypeInd}.spatialRFArray);
    for xcell = 1:nCells(1)
        for ycell = 1:nCells(2)
            hold on; 
            % center
            plot(obj.mosaic{cellTypeInd}.spatialRFcontours{xcell,ycell,1}(1,2:end),...
                obj.mosaic{cellTypeInd}.spatialRFcontours{xcell,ycell,1}(2,2:end),...
                'color',cmap(cellTypeInd,:));
            hold on;
            % surround
            plot(obj.mosaic{cellTypeInd}.spatialRFcontours{xcell,ycell,2}(1,2:end),...
                obj.mosaic{cellTypeInd}.spatialRFcontours{xcell,ycell,2}(2,2:end),...
                'color',cmap(cellTypeInd+8,:));
        end
    end
    axis equal
    title(sprintf('%s',obj.mosaic{cellTypeInd}.nameCellType),'fontsize',16);
    xlabel(sprintf('Distance (\\mum)'),'fontsize',16);
    ylabel(sprintf('Distance (\\mum)'),'fontsize',16);
end

