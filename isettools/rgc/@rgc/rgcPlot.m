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

%% Plot the membrane voltages for a random trial

figure;
szSpike = size(horzcat(obj.mosaic{1}.spikeResponse{1,1,:,2}));
maxTrials = szSpike(2);
for cellTypeInd = 1:obj.numberCellTypes
    subplot(3,2,cellTypeInd);
    plot(horzcat(obj.mosaic{cellTypeInd}.spikeResponse{:,:,floor(maxTrials*rand)+1,2}))
    xlabel(sprintf('Time (msec)'),'fontsize',16);
    ylabel(sprintf('Membrane Voltage (\\muV)'),'fontsize',16);
    title(sprintf('%s',obj.mosaic{cellTypeInd}.nameCellType),'fontsize',16);
end



%% Plot the spike rasters and/or PSTH for each cell

dt = .01; % make this a get from sensor

for cellTypeInd = 1:obj.numberCellTypes
    figure;
    cellCtr = 0;
    clear psth tsp mtsp
    
    nCells = size(obj.mosaic{cellTypeInd}.spatialRFArray);
    
    for xcell = 1:nCells(1)
        for ycell = 1:nCells(2)
            
            cellCtr = cellCtr+1;
            
            for trial = 1:maxTrials
                tsp{trial} = obj.mosaic{cellTypeInd}.spikeResponse{xcell,ycell,trial,1};
            end
            
            subplot(nCells(1),nCells(2),cellCtr);
            % subplot(2*nCells(1),nCells(2),nCells(1)+nCells(1)*(2*(xcell-1))+ycell);
            mtsp = plotraster(tsp);

            [psth{cellCtr},tt,pstv,spr] = compPSTH(mtsp*dt, .001, .002, [0 1], .005);
            axis([0 30 0 maxTrials]);
            
            
            % subplot(nCells(1),nCells(2),cellCtr);
            % plot(tt/.01,psth{cellCtr});
            
            % hold on;
            % gf = obj.mosaic{1}.generatorFunction;
            % plot((obj.mosaic{cellTypeInd}.linearResponse{xcell,ycell}));            
            % axis([0 30 min(obj.mosaic{cellTypeInd}.linearResponse{xcell,ycell}) max(obj.mosaic{cellTypeInd}.linearResponse{xcell,ycell})]);
        end
    end
    
    suptitle(sprintf('%s',obj.mosaic{cellTypeInd}.nameCellType));
end

%% Plot the raster/PSTH from a single cell
szSpike = size(horzcat(obj.mosaic{1}.spikeResponse{1,1,:,2}));
maxTrials = szSpike(2); dt = .01;
cellTypeInd = 4;
nCells = size(obj.mosaic{cellTypeInd}.spatialRFArray);

figure;

xcell = 1; ycell = 1;

for trial = 1:maxTrials
    tsp{trial} = obj.mosaic{cellTypeInd}.spikeResponse{xcell,ycell,trial,1};
end

subplot(2,1,1);
mtsp = plotraster(tsp);
axis([0 30 0 maxTrials]);

xlabel(sprintf('Time (msec)'),'fontsize',16);
ylabel(sprintf('Spikes'),'fontsize',16);
title(sprintf('Raster, cell at x = %d, y = %d, %s',xcell,ycell,obj.mosaic{cellTypeInd}.nameCellType),'fontsize',16);

subplot(2,1,2);
[psth,tt,pstv,spr] = compPSTH(mtsp*dt, .001, .002, [0 1], .005);
plot(tt/.01,psth);
axis([0 30 0 max(psth)]);

xlabel(sprintf('Time (msec)'),'fontsize',16);
ylabel(sprintf('Smoothed Spike Rate (spikes/sec)'),'fontsize',16);
title(sprintf('Smoothed Rate, cell at x = %d, y = %d, %s',xcell,ycell,obj.mosaic{cellTypeInd}.nameCellType),'fontsize',16);
     