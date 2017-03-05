function [uData, hf] = plot(obj, plotType, varargin)
% Plot function for rgcMosaic
%    [uData, hf] = rgcMosaic.plot(plotType, varargin)
%
% Required input
%   plotType - string, type of plot
%
% Optional input (key-val pairs in varargin):
%   'hf' - figure or axis handle
%      By default, this is the figure handle of the rgcMosaic window
%      It can also be one of the two axes in the window.
%      Otherwise
%             []            - create plot in new figure
%             figure handle - an alternative figure
%             'none'        - don't plot, only generate uData
%
% Outputs:
%   uData - Computed user data
%   hf    - figure handle
%
% Plot type can be chosen from
%   'spike mean image'    - Gray scale image of mean spikes
%   'mosaic'              - RGC mosaic as circles
%   'linear movie'        - Linear voltages as movie
%   'psth'                - Peristimulus time histogram as ...
%
% Example:
%
% BW, ISETBIO TEAM, 2016

%% parse input
p = inputParser;
p.KeepUnmatched = true;

% What about obj?
p.addRequired('plotType', @isstr);                        % Type of plot
p.addParameter('hf', obj.figureHandle, @isgraphics);  % figure handle

p.parse(plotType, varargin{:});
hf = p.Results.hf;

uData = [];

% Select place for plot to appear
% The mosaic window has a response image and a geometry image
if isempty(hf), hf = vcNewGraphWin;
elseif isgraphics(hf, 'figure'), figure(hf);
elseif isgraphics(hf, 'axes'), axes(hf);
end

fprintf('Plot %s for %s class\n',plotType,class(obj));

switch ieParamFormat(plotType)
    case 'spikemeanimage'
        % Spike mean image
        g = guidata(hf);
        axes(g.axisResponse);
        
        spikes = obj.get('spikes');
        img = mean(spikes,3);
        colormap(gray(256)); imagesc(img); axis image;
        set(gca,'xticklabels','','yticklabels','');
        colorbar; drawnow;
        
    case 'spikemovie'
        % Movie of spiking activity
        % Should this work with the mosaic plot of circles?
        psthTest = obj.get('spikes');
        clear vParams; vParams = [];
        vParams.FrameRate = 30; vParams.show = true; %vParams.step = 2;
        frameSkip = round(1./obj.get('dt'));
        ieMovie(psthTest(:,:,1:frameSkip:end),vParams);
        
    case{'linearmovie'}
        % Continuous voltages prior to spike generation
        responseLinear = obj.get('responseLinear');
        
        clear vParams; 
        vParams.FrameRate = 30; 
        vParams.show = true;
        
        % Should we add controls like the cone mosaic window, or just play
        % it once?
        ieMovie(responseLinear,vParams);
      
        
    case 'psthmeanimage'
        psth = obj.get('psth');
        imagesc(mean(psth,3));
        axis image; colormap(gray(256)); colorbar;
        set(gca,'xticklabels','','yticklabels','');
        drawnow;

    case 'psth'
        % Peri-stimulus time graph for all cells.
        % Kind of a weird plot to make.
        timeStep = obj.dt;
        psth = obj.get('psth');
        resp = RGB2XWFormat(psth);
        
        g = guidata(hf);
        axes(g.axisResponse);
        cla(g.axisResponse,'reset');
        plot(g.axisResponse,timeStep*(1:size(resp,1)),sum(resp,2));
        xlabel('Time (sec)');
        ylabel(sprintf('Spikes per %d ms',timeStep*1e3));
        grid on; axis image
        
    case 'mosaic'
        % Plot the mosaic spatial receptive field geometry
        % We need to add the possibility of elliptical forms some day.
        % And we should figure out how to do center/surround
        cla reset;
        
        % Oddly, the center is (row,col)
        center = cell2mat(obj.cellLocation(:));  % um w.r.t. center of image
        radius = obj.rfDiameter/2;
        ieShape('circle','center',center,...
            'radius',radius, ...
            'color','b');
        
        set(gca,...
            'xlim',[min(center(:,2)) - radius, max(center(:,2)) + radius],...
            'ylim',[min(center(:,1)) - radius, max(center(:,1)) + radius]);
        xlabel(sprintf('Distance (\\mum)'),'fontsize',14);
        
    otherwise
        error('Unknown plot type %s\n',plotType);
end

end

%%
% Plot type can be chosen from
%         'rf',...              - (center - surround) spatial RF surfaces
%         'rfImage',...         - a (center - surround) spatial RF image
%         'mosaic',...          - the 1 STD spatial RF mosaic of each type
%         'sRFcenter',...       - center spatial RF surfaces
%         'sRFsurround',...     - surround spatial RF surfaces
%         'temporal',...        - (center - surround) temporal impulse responses
%         'tCenter',...         - center temporal impulse response
%         'tSurround',...       - surround temopral impulse response
%         'postSpikeFilter',... - post-spike filter time course
%         'couplingFilter',...  - coupling filters time course
%         'responseLinear',...  - linear response of all cells
%         'nlResponse',...      - nonlinear response fgenerator(linear) of all cells
%         'responseSpikes',...   - average waveform over N trials including
%                                   post-spike and coupling filter effects
%         'raster',...          - spike rasters of all cells from N trials
%         'psth'...             - peristimulus time histogram (specify a cell)
%
% Examples:
%   irPlot(innerRetina,'mosaic');
%   irPlot(innerRetina,'psth');
%   irPlot(innerRetina,'psth','type','onParasol');
%   irPlot(innerRetina,'psth','cell',[1 1]);
%   irPlot(innerRetina,'psth','type','onParasol','cell',[1 1]);
%
% JRG/HJ/BW, ISETBIO TEAM, 2016
%
%
% % parse input
% p = inputParser;
% p.KeepUnmatched = true;
% p.addRequired('type', @isstr);               % Type of plot
% p.addParameter('hf', []);                    % figure handle
% p.addParameter('oi',[],@isstruct);           % Used for spectral qe
% p.addParameter('mosaic',[1:length(obj.mosaic)],@isnumeric); % specify mosaic
% p.addParameter('cell',[],@isnumeric);
% p.addParameter('dt',1,@isnumeric);
% p.addParameter('color','b',@ischar);
% p.parse(type, varargin{:});
% hf = p.Results.hf;
% oi = p.Results.oi;
% mosaicIndices = p.Results.mosaic;
% cellIndices = p.Results.cell;
% dt = p.Results.dt;
% color = p.Results.color;
% uData = [];
%
% % plot
% if isempty(hf), hf = vcNewGraphWin;
% elseif isgraphics(hf, 'figure'), figure(hf);
% elseif isgraphics(hf, 'axes'), axes(hf);
% end
%
% switch ieParamFormat(type)
%
%     case{'ecc'}
%         %%% The temporal equivalent eccentricity of the retinal patch
%         plotPatchEccentricity(obj.eyeAngle, obj.eyeRadius, obj.eyeSide, obj.temporalEquivEcc);
%
%     case{'mosaic'}
%         %%% Plot the mosaic of each RGC type
%
%         for cellTypeInd = mosaicIndices
%
%             % Get contour lines for mosaic RFs
%             spatialRFcontours = plotContours(obj.mosaic{cellTypeInd}, obj.spacing, obj.col);
%
%             % Subplot if more than one mosaic is being plotted
%             if length(mosaicIndices)>1; subplot(ceil(length(mosaicIndices)/2),2,cellTypeInd); end;
%
%             nCells = obj.mosaic{cellTypeInd}.get('mosaic size');
%
%             % Convert RGC position to distance
%             patchSizeX = obj.spacing;
%             numberCellsX = obj.col;
%             umPerCell = 1e6*patchSizeX/numberCellsX;
%
%             cmap = parula(16);
%             for xcell = 1:nCells(1)
%                 for ycell = 1:nCells(2)
%                     hold on;
%                     % center
%                     plot(umPerCell*spatialRFcontours{xcell,ycell,1}(1,2:end),...
%                         umPerCell*spatialRFcontours{xcell,ycell,1}(2,2:end),...
%                         'color',cmap(cellTypeInd,:));
%                     hold on;
%                     % surround
%                     plot(umPerCell*spatialRFcontours{xcell,ycell,2}(1,2:end),...
%                         umPerCell*spatialRFcontours{xcell,ycell,2}(2,2:end),...
%                         'color',cmap(cellTypeInd+8,:));
%                 end
%             end
%             axis equal
%             title(sprintf('%s',obj.mosaic{cellTypeInd}.cellType),'fontsize',14);
%             xlabel(sprintf('Distance (\\mum)'),'fontsize',14);
%             ylabel(sprintf('Distance (\\mum)'),'fontsize',14);
%
%         end
%
%     case{'rf'}
%         %%% A surface representing the RF (center - surround)
%         for cellTypeInd = mosaicIndices
%
%              % Subplot if more than one mosaic is being plotted
%             if length(obj.mosaic)>1; subplot(ceil(length(obj.mosaic)/2),2,cellTypeInd); end;
%
%             % Convert RGC position to distance
%             patchSizeX = obj.spacing;
%             numberCellsX = obj.col;
%             umPerCell = 1e6*patchSizeX/numberCellsX;
%
%             % Get x and y axes in distance
%             xd = umPerCell*([1:size(obj.mosaic{cellTypeInd}.sRFcenter{1,1},1)]);
%             yd = umPerCell*([1:size(obj.mosaic{cellTypeInd}.sRFcenter{1,1},2)]);
%             % Take difference of center and surround
%             rfSurf = obj.mosaic{cellTypeInd}.sRFcenter{1,1}-obj.mosaic{cellTypeInd}.sRFsurround{1,1};
%             surface(xd,yd,rfSurf); shading flat; view(40,40);
%
%             title(sprintf('Spatial Receptive Field, %s',obj.mosaic{cellTypeInd}.cellType),'fontsize',14);
%             xlabel(sprintf('Distance (\\mum)'),'fontsize',14);
%             ylabel(sprintf('Distance (\\mum)'),'fontsize',14);
%             zlabel(sprintf('Response (spikes/sec)'),'fontsize',14);
%
%             surfMin = -max(obj.mosaic{cellTypeInd}.sRFsurround{1,1}(:));
%             surfMax = max(obj.mosaic{cellTypeInd}.sRFcenter{1,1}(:));
%             axis([0 xd 0 yd  surfMin surfMax ]);
%         end
%
%     case{'temporal','temporalfilter'}
%         %%% Plot the RGB impulse response of each mosaic
%         for cellTypeInd = mosaicIndices
%
%             % Subplot if more than one mosaic is being plotted
%             if length(mosaicIndices)>1; subplot(ceil(length(mosaicIndices)/2),2,cellTypeInd); end;
%
%             % Get time and filter data
%             tFilter = (obj.mosaic{cellTypeInd}.tCenter{1,1});
%             tx = obj.timing:obj.timing:obj.timing*length(tFilter);
%             if length(tFilter)>1
%                 plot(tx,tFilter);
%
%             else
%                 scatter(tx,tFilter,'x');
%                 disp('Temporal filter is identity, see bipolar temporal filter');
%             end
%
%             title(sprintf('Temporal Filter, RGB, %s',obj.mosaic{cellTypeInd}.cellType),'fontsize',14);
%             xlabel(sprintf('Time (sec)'),'fontsize',14);
%             ylabel(sprintf('Conditional Intensity (AU)'),'fontsize',14);
%         end
%
%     case{'couplingfilter','coupling'}
%
%         %%% Plot the post spike filter of each cell type
%         for cellTypeInd = mosaicIndices
%
%             % Subplot if more than one mosaic is being plotted
%             if length(mosaicIndices)>1; subplot(ceil(length(mosaicIndices)/2),2,cellTypeInd); end;
%
%             if iscell(obj.mosaic{cellTypeInd}.couplingFilter)
%                 cplf = (RGB2XWFormat(obj.mosaic{cellTypeInd}.couplingFilter{1,1}));
%                 tx = obj.timing:obj.timing:obj.timing*length(cplf);
%                 % Subplot if more than one mosaic is being plotted
%                 if length(mosaicIndices)>1; subplot(ceil(length(mosaicIndices)/2),2,cellTypeInd); end;
%
%                 plot(tx, (cplf'));
%                 title(sprintf('Coupling Filters, %s',obj.mosaic{cellTypeInd}.cellType),'fontsize',14);
%                 xlabel(sprintf('Time (sec)'),'fontsize',14);
%                 ylabel(sprintf('Response (spikes/sec)'),'fontsize',14);
%
%             end
%         end
%
%     case{'linear', 'responselinear'}
%         timeStep = obj.timing;
%
%         for cellTypeInd = mosaicIndices
%
%             nCells = obj.mosaic{cellTypeInd}.get('mosaic size');
%
%             % Subplot if more than one mosaic is being plotted
%             if length(mosaicIndices)>1; subplot(ceil(length(mosaicIndices)/2),2,cellTypeInd); end;
%
%             resp = RGB2XWFormat(obj.mosaic{cellTypeInd}.responseLinear);
%
%             plot(timeStep:timeStep:timeStep*size(resp,2),resp');
%             xlabel(sprintf('Time (sec)'),'fontsize',14);
%             ylabel(sprintf('Conditional Intensity'),'fontsize',14);
%             title(sprintf('%s',obj.mosaic{cellTypeInd}.cellType),'fontsize',14);
%
%         end
%     case{'movielinear','linearmovie'}
%
%         resp = obj.mosaic{mosaicIndex}.get('response linear');
%         ieMovie(resp);
%
%     case{'raster','rasterresponse'}
%         % Plot spike rasterse
%
%         timeStep = obj.timing*dt;
%         % bindur = ir.mosaic{1}.get('dt');
%
%         for cellTypeInd = mosaicIndices
%
%             nCells = obj.mosaic{cellTypeInd}.get('mosaic size');
%             nCellStart = [1 1]; nCellEnd = nCells;
%             if length(cellIndices)>0
%                 nCellStart = [cellIndices];
%                 nCellEnd = [cellIndices];
%             end
%             nTrials = obj.mosaic{1}.get('numbertrials');
%
%             lastspiketime = obj.mosaic{1}.get('lastspiketime');
%             cellCtr = 0;
%
%             for xcell = nCellStart(1):nCellEnd(1)
%                 for ycell = nCellStart(2):nCellEnd(2)
%                     cellCtr = cellCtr+1;
%                     if length(cellIndices)>2
%                         subplot(length(nCellStart(1):nCellEnd(1)),length(nCellStart(2):nCellEnd(2)),cellCtr);
%                     end
%                     for tr = 1:nTrials
%                         % Get spike times
%                         spikeTimes = obj.mosaic{cellTypeInd}.responseSpikes{xcell,ycell,tr};
%                         hold on;
%
%                         % Draw raster plots
%                         % line([spikeTimes, spikeTimes]*timeStep,[tr tr-1],'color','k');
%                         scatter([spikeTimes].*timeStep,[tr*ones(length(spikeTimes),1)],8,'o',color,'filled');
%
%                         axis([0 timeStep*lastspiketime 0 nTrials+1]);
%
%                     end
%
%                     if cellCtr == 1
%
%                         xlabel(sprintf('Time (sec)'));
%                         ylabel(sprintf('Trial'));
%                         title(sprintf('%s cell [%d %d] raster',obj.mosaic{cellTypeInd}.cellType,xcell,ycell));
%                     end
%
%                 end
%             end
%         end
%
%     case{'psth'}
%         timeStep = obj.timing;
%
%         for cellTypeInd = mosaicIndices
%
%         nCells = obj.mosaic{cellTypeInd}.get('mosaic size');
%             nCellStart = [1 1]; nCellEnd = nCells;
%             if length(cellIndices)>0
%                 nCellStart = [cellIndices];
%                 nCellEnd   = [cellIndices];
%             end
%             nTrials = obj.mosaic{1}.get('numbertrials');
%
%             lastspiketime = obj.mosaic{1}.get('lastspiketime');
%             cellCtr = 0;
%
%             for xcell = nCellStart(1):nCellEnd(1)
%                 for ycell = nCellStart(2):nCellEnd(2)
%                     cellCtr = cellCtr+1;
%
%                     psth = obj.mosaic{cellTypeInd}.get('psth');
%
%                     resp = RGB2XWFormat(psth(xcell,ycell,:));
%
%                     plot(timeStep:timeStep:timeStep*size(resp,2),resp);
%
%
%                     if cellCtr == 1
%                         xlabel(sprintf('Time (sec)'));
%                         ylabel(sprintf('Conditional Intensity'));
%                         title(sprintf('%s',obj.mosaic{cellTypeInd}.cellType));
%                     end
%                 end
%             end
%         end
% end
%
