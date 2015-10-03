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
plotPatchEccentricity(obj.patchLocationPolarAngleDegrees, obj.patchLocationPolarRadiusMicrometers, obj.eyeLeftOrRight, obj.temporalEquivEcc);

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

%% Plot the spatial and temporal impulse responses

% cellTypeInd = 1;

% figure; plot(.01:.01:.2,bsxfun(@plus,horzcat(rgc1.mosaic{cellTypeInd}.temporalImpulseResponseCenterRGB{:}),[0 0 0.01]))
% title(sprintf('Temporal Impulse Response, RGB, %s',rgc1.mosaic{cellTypeInd}.nameCellType),'fontsize',16);
% xlabel(sprintf('Time (sec)'),'fontsize',16);
% ylabel(sprintf('Response (spikes/sec)'),'fontsize',16);
% 
% psf = squeeze(rgc1.mosaic{cellTypeInd}.couplingFilter{1,1}(1,1,:));
% figure; plot((1:length(psf))./1000, psf);
% title(sprintf('Exponentiated Post-Spike Filter, %s',rgc1.mosaic{cellTypeInd}.nameCellType),'fontsize',16);
% xlabel(sprintf('Time (sec)'),'fontsize',16);
% ylabel(sprintf('Response (spikes/sec)'),'fontsize',16);
% 
% 
% cplf = (horzcat(rgc1.mosaic{cellTypeInd}.couplingFilter{:}));
% figure; plot((1:length(cplf))./1000, squeeze(cplf(1,:,:)));
% title(sprintf('Exponentiated Coupling Filters, %s',rgc1.mosaic{cellTypeInd}.nameCellType),'fontsize',16);
% xlabel(sprintf('Time (sec)'),'fontsize',16);
% ylabel(sprintf('Response (spikes/sec)'),'fontsize',16);
% 
% figure; surface(rgc1.mosaic{cellTypeInd}.spatialRFArray{1,1}); shading flat;
% title(sprintf('Spatial Receptive Field, %s',rgc1.mosaic{cellTypeInd}.nameCellType),'fontsize',16);
% xlabel(sprintf('Distance (\\mum)'),'fontsize',16);
% ylabel(sprintf('Distance (\\mum)'),'fontsize',16);
% zlabel(sprintf('Response (spikes/sec)'),'fontsize',16);
%% Plot the membrane voltages for a random trial

figure;
szSpike = size(horzcat(obj.mosaic{3}.spikeResponse{1,1,:,2}));
maxTrials = szSpike(2);
for cellTypeInd = 3%1:obj.numberCellTypes
%     subplot(3,2,cellTypeInd);
    plot(horzcat(obj.mosaic{cellTypeInd}.spikeResponse{:,:,floor(maxTrials*rand)+1,2}))
    xlabel(sprintf('Time (msec)'),'fontsize',16);
    ylabel(sprintf('Membrane Voltage (\\muV)'),'fontsize',16);
    title(sprintf('%s',obj.mosaic{cellTypeInd}.nameCellType),'fontsize',16);
end



%% Plot the spike rasters and/or PSTH for each cell

dt = .01; % make this a get from sensor

for cellTypeInd = 3%1:obj.numberCellTypes
    figure;
    cellCtr = 0;
    clear psth tsp mtsp
    
    nCells = size(obj.mosaic{cellTypeInd}.spatialRFArray);
    
    rasterResponse =  mosaicGet(obj.mosaic{cellTypeInd}, 'rasterResponse');
    
    for xcell = 1:nCells(1)
        for ycell = 1:nCells(2)
            
            cellCtr = cellCtr+1;
            
%             if ~sum(cellfun(@isempty,tsp))
            
            for trial = 1:maxTrials
                tsp{trial} = obj.mosaic{cellTypeInd}.spikeResponse{xcell,ycell,trial,1};
            end
            
            subplot(nCells(1),nCells(2),cellCtr);
            % subplot(2*nCells(1),nCells(2),nCells(1)+nCells(1)*(2*(xcell-1))+ycell);
            
            if sum(cellfun(@isempty,tsp))~=maxTrials
            
                mtsp = plotraster(tsp);
            else
                mtsp = [];
            end
            raster{xcell,ycell} = mtsp;
                        
            % [psth{cellCtr},tt,pstv,spr] = compPSTH(mtsp*dt, .001, .002, [0 1], .005);
            [psth{xcell,ycell},tt,pstv,spr] = compPSTH(mtsp*dt, .001, .002, [0 1], .005);
            axis([0 30 0 maxTrials]);
            
            
            % subplot(nCells(1),nCells(2),cellCtr);
%             plot(tt/.01,psth{xcell,ycell});
%             if ~isnan(psth{xcell,ycell})
%                 axis([0 30 0 max(psth{xcell,ycell})]);
%             end

            % hold on;
            % gf = obj.mosaic{1}.generatorFunction;
            % plot((obj.mosaic{cellTypeInd}.linearResponse{xcell,ycell}));            
            % axis([0 30 min(obj.mosaic{cellTypeInd}.linearResponse{xcell,ycell}) max(obj.mosaic{cellTypeInd}.linearResponse{xcell,ycell})]);
            
        end
    end
    % mosaicSet(obj.mosaic{cellTypeInd},'rasterResponse',raster);
    % mosaicSet(obj.mosaic{cellTypeInd},'psthResponse',psth);
    suptitle(sprintf('%s',obj.mosaic{cellTypeInd}.nameCellType));
end

%% Plot the raster/PSTH from a single cell
% szSpike = size(horzcat(obj.mosaic{1}.spikeResponse{1,1,:,2}));
% maxTrials = szSpike(2); dt = .01;
% cellTypeInd = 4;
% nCells = size(obj.mosaic{cellTypeInd}.spatialRFArray);

figure;

xcell = 4; ycell = 2;
cellTypeInd = 3;
for trial = 1:maxTrials
    tsp{trial} = obj.mosaic{cellTypeInd}.spikeResponse{xcell,ycell,trial,1};
end

subplot(2,1,1);
            
if sum(cellfun(@isempty,tsp))~=maxTrials
    
    mtsp = plotraster(tsp);
else
    mtsp = [];
end

axis([0 30 0 maxTrials]);

xlabel(sprintf('Time (msec)'),'fontsize',16);
ylabel(sprintf('Spikes'),'fontsize',16);
title(sprintf('Raster, cell at x = %d, y = %d, %s',xcell,ycell,obj.mosaic{cellTypeInd}.nameCellType),'fontsize',16);

subplot(2,1,2);
[psth,tt,pstv,spr] = compPSTH(mtsp*dt, .001, .002, [0 1], .005);
plot(tt/.01,psth);
if ~isnan(psth)
    axis([0 30 0 .0001+max(psth)]);
end

xlabel(sprintf('Time (msec)'),'fontsize',16);
ylabel(sprintf('Smoothed Spike Rate (spikes/sec)'),'fontsize',16);
title(sprintf('Smoothed Rate, cell at x = %d, y = %d, %s',xcell,ycell,obj.mosaic{cellTypeInd}.nameCellType),'fontsize',16);
     