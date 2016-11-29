function [uData, hf] = plot(obj, type, varargin)
% Plot function for rgcMosaic
%    [uData, hf] = rgcMosaic.plot(type, varargin)
%
% Required input
%   type - string, type of plot
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
%   'rgc mosaic'          - Color image of the cone arrangement
%   'spike mean image'     - Gray scale movie of current
%   'spike movie'   - Cone photocurrent graphs
%
% Example:
%
% BW, ISETBIO TEAM, 2016

%% parse input
p = inputParser;
p.KeepUnmatched = true;

% What about obj?
p.addRequired('type', @isstr);                        % Type of plot
p.addParameter('hf', obj.figureHandle, @isgraphics);  % figure handle

p.parse(type, varargin{:});
hf = p.Results.hf;

uData = [];

% Select place for plot to appear
% The mosaic window has a response image and a geometry image
if isempty(hf), hf = vcNewGraphWin;  
elseif isgraphics(hf, 'figure'), figure(hf); 
elseif isgraphics(hf, 'axes'), axes(hf);
end

% % set color order so that LMS plots as RGB
% if ~isequal(hf, 'none')
%     co = get(gca, 'ColorOrder');
%     if isgraphics(hf,'axes')
%         set(get(hf,'parent'),'DefaultAxesColorOrder', co([2 5 1 3 4 6 7], :)) 
%     else  % Figure
%         set(hf, 'DefaultAxesColorOrder', co([2 5 1 3 4 6 7], :));
%     end
% end

fprintf('Plot %s for %s class\n',type,class(obj));

switch ieParamFormat(type)
    case 'spikemeanimage'
        % Spike mean image
        g = guidata(hf);
        axes(g.axisResponse);
        
        spikes = obj.get('responseSpikes');
        img = mean(spikes,3);
        imagesc(img);
        
    case 'spikemovie'
        disp('Spike Movie')
    case{'movielinear','linearmovie'}
        resp = obj.get('response linear');
        dFlag.FrameRate = 2;
        g = guidata(hf);
        axes(g.axisResponse);
        
        ieMovie(resp,'dFlag',dFlag);
    case 'psth'
        % Peri-stimulus time histogram
        timeStep = obj.Parent.timing;
        psth = obj.get('psth');
        resp = RGB2XWFormat(psth);
        
        g = guidata(hf);
        axes(g.axisResponse); 
        cla(g.axisResponse,'reset');
        plot(g.axisResponse,timeStep:timeStep:timeStep*size(resp,1),sum(resp,2));
        xlabel(sprintf('Time (sec)'));
        % ylabel(sprintf('Conditional Intensity'));
        % title(sprintf('%s',obj.cellType));
        
        
    case 'mosaic'
        % Plot the mosaic spatial receptive field geometry
        % We should plot, say, 5x5 array just to illustrate, rather than
        % the whole mosaic.  Only the whole mosaic on demand
        g = guidata(obj.figureHandle);
        cla reset;
        
        % Somehow, we need these variables, too.  Let's rethink the
        % parameterization of the spatial receptive fields.
        %
        %         metersPerPixel = (spacing/1e-6)/col;
        %         rfPixels = obj.rfDiameter/metersPerPixel;
        %         extent = .5*round(size(obj.sRFcenter{1,1},1)/rfPixels);
        %

        % Get contour lines for mosaic RFs
        % See below:  Replace this with ieShape()
        % spatialRFcontours = plotContours(obj, obj.Parent.spacing, obj.Parent.col);
        
        
        % Convert RGC position to distance
        % These calculations can't be right.  The RF size has to be
        % independent of the number of cells.... FIX.
        patchSizeX      = obj.Parent.spacing;   % Meters
        numberBipolarsX = obj.Parent.col;       % Bipolar cells
        umPerBipolar       = 1e6*patchSizeX/numberBipolarsX;   % Meters per Bipolar
        
        % Replace with ieShape('circle','center',...,'radius',...)
        % To draw the center and surround mosaic geometry
        % Axis units should be um
        
        % New strategy: make surface plot of all RFs combined and take a
        % thin z slice to show contours. This is useful for when RFs are
        % not circular, and each has a different elliptical shape (JRG).
        
        nCells = obj.get('mosaic size');
%         nCell = 5;   % Central 5 cells
%         cellList = -2:2;
        % I am noticing that the center of the sRF when the size is even
        % is, well, not good.  We should probably always insist on odd
        % row/col of the sRF images.
        for xcell = 1:nCells(1)
            for ycell = 1:nCells(2)
%         for xcell = cellList
%             for ycell = cellList
%                 center = [xcell*umPerCell, ycell*umPerCell];
%                 % This should be a real radius ... it appears to be just
%                 % the spacing ... must read the constructor (BW).
%                 radius = umPerCell;   
%                 % Could be out of the loop ... or could be pulled from a
%                 % space varying rgc object.
%                 [h,pts] = ieShape('circle','center',center,'radius',radius,'color','b');
%                 fill(pts(:,1),pts(:,2),[1 0 1]);
%                 % center
%                 %                 plot(umPerCell*spatialRFcontours{xcell,ycell,1}(1,2:end),...
%                 %                     umPerCell*spatialRFcontours{xcell,ycell,1}(2,2:end),...
%                 %                     'color','b');
%                 %                 hold on;
%                 %                 % surround
%                 %                 plot(umPerCell*spatialRFcontours{xcell,ycell,2}(1,2:end),...
%                 %                     umPerCell*spatialRFcontours{xcell,ycell,2}(2,2:end),...
%                 %                     'color','y');
                
                hold on;
                % Generate x and y coordinates for spatial RF
                % These are in units of bipolars ["the RGC RF is N bipolar
                % wide"], converted to meters below
                sRFx = [1:size(obj.sRFcenter{xcell,ycell},2)];
                sRFy = [1:size(obj.sRFcenter{xcell,ycell},1)];
                
                % Shift to zero by subtracting half of size
                sRFxZero = round(size(obj.sRFcenter{xcell,ycell},2)/2);
                sRFyZero = round(size(obj.sRFcenter{xcell,ycell},1)/2);
                
                % Get center coordinate
                sRFxCenter = obj.cellLocation{xcell,ycell}(2);
                sRFyCenter = obj.cellLocation{xcell,ycell}(1);
                
                % Combine for appropriate coordinates
                plotX = sRFx-sRFxZero+sRFxCenter;
                plotY = sRFy-sRFyZero+sRFyCenter;
                
                % Plot surface
                % Convert coordinates from number of bipolar cells to
                % meters
                surf(umPerBipolar*plotX,umPerBipolar*plotY,obj.sRFcenter{xcell,ycell});
        
            end
        end
        
        axis equal
        shading interp
        
        % Since we've plotted the surface, take a thin z slice for contours
        maxRF = max(obj.sRFcenter{1,1}(:));
        ax1 = axis
        axis([ax1 exp(-1)*2*maxRF-.05 exp(-1)*2*maxRF]);
%         axis([-10-max(umPerBipolar*plotX) 10+max(umPerBipolar*plotX) -10-max(umPerBipolar*plotY) 10+max(umPerBipolar*plotY) .5*maxRF-.05 .5*maxRF]);
        
        
%         alim = 1.5*[-nCell*umPerCell, nCell*umPerCell]/2;
%         set(gca,'xlim',alim,'ylim',alim);
        xlabel(sprintf('Distance (\\mum)'),'fontsize',14);
%         grid on
        hold off;
    otherwise
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
