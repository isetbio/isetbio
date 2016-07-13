function [hf, uData] = irPlot(obj, type, varargin)
% Plot function for ir
%    [hf, uData] = ir.plot(type, varargin)
%
% Inputs:
%   type - string, type of plot
%
% Optional input (key-val pairs in varargin):
%   'hf' - figure handle or control structure, the meaning of value is
%             []            - create plot in new figure
%             figure handle - plot in figure specified by hf
%             'none'        - don't plot
%
% Outputs:
%   hf    - figure handle
%   uData - Computed data
%
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
%         'rasterResponse',...  - spike rasters of all cells from N trials
%         'psthResponse'...     - peristimulus time histogram responses of all cells
%
% Examples:
%   irPlot(innerRetina,'mosaic');
%   irPlot(innerRetina,'psth');
%   irPlot(innerRetina,'psth','type','onParasol');
%   irPlot(innerRetina,'psth','cell',[1 1]);
%   irPlot(innerRetina,'psth','type','onParasol','cell',[1 1]);
%
% JRG/HJ/BW, ISETBIO TEAM, 2016


% parse input
p = inputParser;
p.KeepUnmatched = true;
p.addRequired('type', @isstr);               % Type of plot
p.addParameter('hf', []);                    % figure handle
p.addParameter('oi',[],@isstruct);           % Used for spectral qe
p.addParameter('mosaic',[1:length(obj.mosaic)],@isnumeric); % specify mosaic 
p.addParameter('cell',[],@isnumeric);
p.addParameter('dt',1,@isnumeric);
p.addParameter('color','b',@ischar);
p.parse(type, varargin{:});
hf = p.Results.hf;
oi = p.Results.oi;
mosaicIndices = p.Results.mosaic;
cellIndices = p.Results.cell;
dt = p.Results.dt;
color = p.Results.color;
uData = [];

% plot
if isempty(hf), hf = vcNewGraphWin;
elseif isgraphics(hf, 'figure'), figure(hf);
elseif isgraphics(hf, 'axes'), axes(hf);
end

switch ieParamFormat(type)   
        
    case{'ecc'}
        %%% The temporal equivalent eccentricity of the retinal patch
        plotPatchEccentricity(obj.eyeAngle, obj.eyeRadius, obj.eyeSide, obj.temporalEquivEcc);
        
    case{'mosaic'}
        %%% Plot the mosaic of each RGC type
        
        for cellTypeInd = mosaicIndices
            
            % Get contour lines for mosaic RFs
            spatialRFcontours = plotContours(obj.mosaic{cellTypeInd}, obj.spacing, obj.col);
            
            % Subplot if more than one mosaic is being plotted
            if length(mosaicIndices)>1; subplot(ceil(length(mosaicIndices)/2),2,cellTypeInd); end;
            
            nCells = obj.mosaic{cellTypeInd}.get('mosaic size');
            
            % Convert RGC position to distance
            patchSizeX = obj.spacing;
            numberCellsX = obj.col;            
            umPerCell = 1e6*patchSizeX/numberCellsX;
            
            cmap = parula(16);
            for xcell = 1:nCells(1)
                for ycell = 1:nCells(2)
                    hold on;
                    % center
                    plot(umPerCell*spatialRFcontours{xcell,ycell,1}(1,2:end),...
                        umPerCell*spatialRFcontours{xcell,ycell,1}(2,2:end),...
                        'color',cmap(cellTypeInd,:));
                    hold on;
                    % surround
                    plot(umPerCell*spatialRFcontours{xcell,ycell,2}(1,2:end),...
                        umPerCell*spatialRFcontours{xcell,ycell,2}(2,2:end),...
                        'color',cmap(cellTypeInd+8,:));
                end
            end
            axis equal
            title(sprintf('%s',obj.mosaic{cellTypeInd}.cellType),'fontsize',14);
            xlabel(sprintf('Distance (\\mum)'),'fontsize',14);
            ylabel(sprintf('Distance (\\mum)'),'fontsize',14);
          
        end
        
    case{'rf'}
        %%% A surface representing the RF (center - surround)
        for cellTypeInd = mosaicIndices
            
             % Subplot if more than one mosaic is being plotted
            if length(obj.mosaic)>1; subplot(ceil(length(obj.mosaic)/2),2,cellTypeInd); end;
            
            % Convert RGC position to distance
            patchSizeX = obj.spacing;
            numberCellsX = obj.col;
            umPerCell = 1e6*patchSizeX/numberCellsX;
            
            % Get x and y axes in distance
            xd = umPerCell*([1:size(obj.mosaic{cellTypeInd}.sRFcenter{1,1},1)]);
            yd = umPerCell*([1:size(obj.mosaic{cellTypeInd}.sRFcenter{1,1},2)]);
            % Take difference of center and surround
            rfSurf = obj.mosaic{cellTypeInd}.sRFcenter{1,1}-obj.mosaic{cellTypeInd}.sRFsurround{1,1};
            surface(xd,yd,rfSurf); shading flat; view(40,40);
            
            title(sprintf('Spatial Receptive Field, %s',obj.mosaic{cellTypeInd}.cellType),'fontsize',14);
            xlabel(sprintf('Distance (\\mum)'),'fontsize',14);
            ylabel(sprintf('Distance (\\mum)'),'fontsize',14);
            zlabel(sprintf('Response (spikes/sec)'),'fontsize',14);
            
            surfMin = -max(obj.mosaic{cellTypeInd}.sRFsurround{1,1}(:));
            surfMax = max(obj.mosaic{cellTypeInd}.sRFcenter{1,1}(:));
            axis([0 xd 0 yd  surfMin surfMax ]);
        end
        
    case{'temporal','temporalfilter'}
        %%% Plot the RGB impulse response of each mosaic
        for cellTypeInd = mosaicIndices
            
            % Subplot if more than one mosaic is being plotted
            if length(mosaicIndices)>1; subplot(ceil(length(mosaicIndices)/2),2,cellTypeInd); end;           

            % Get time and filter data
            tFilter = (obj.mosaic{cellTypeInd}.tCenter{1,1});
            tx = obj.timing:obj.timing:obj.timing*length(tFilter);
            if length(tFilter)>1
                plot(tx,tFilter);

            else
                scatter(tx,tFilter,'x');
                disp('Temporal filter is identity, see bipolar temporal filter');
            end
            
            title(sprintf('Temporal Filter, RGB, %s',obj.mosaic{cellTypeInd}.cellType),'fontsize',14);
            xlabel(sprintf('Time (sec)'),'fontsize',14);
            ylabel(sprintf('Conditional Intensity (AU)'),'fontsize',14);
        end
        
    case{'couplingfilter','coupling'}
        
        %%% Plot the post spike filter of each cell type
        for cellTypeInd = mosaicIndices
            
            % Subplot if more than one mosaic is being plotted
            if length(mosaicIndices)>1; subplot(ceil(length(mosaicIndices)/2),2,cellTypeInd); end;           

            if iscell(obj.mosaic{cellTypeInd}.couplingFilter)
                cplf = (RGB2XWFormat(obj.mosaic{cellTypeInd}.couplingFilter{1,1}));
                tx = obj.timing:obj.timing:obj.timing*length(cplf);
                % Subplot if more than one mosaic is being plotted
                if length(mosaicIndices)>1; subplot(ceil(length(mosaicIndices)/2),2,cellTypeInd); end;
                
                plot(tx, (cplf'));
                title(sprintf('Coupling Filters, %s',obj.mosaic{cellTypeInd}.cellType),'fontsize',14);
                xlabel(sprintf('Time (sec)'),'fontsize',14);
                ylabel(sprintf('Response (spikes/sec)'),'fontsize',14);
                
            end
        end

    case{'linear', 'responselinear'}  
        timeStep = obj.timing;
                
        for cellTypeInd = mosaicIndices
            
            nCells = obj.mosaic{cellTypeInd}.get('mosaic size');
            
            % Subplot if more than one mosaic is being plotted
            if length(mosaicIndices)>1; subplot(ceil(length(mosaicIndices)/2),2,cellTypeInd); end;           

            resp = RGB2XWFormat(obj.mosaic{cellTypeInd}.responseLinear);
            
            plot(timeStep:timeStep:timeStep*size(resp,2),resp');
            xlabel(sprintf('Time (sec)'),'fontsize',14);
            ylabel(sprintf('Conditional Intensity'),'fontsize',14);
            title(sprintf('%s',obj.mosaic{cellTypeInd}.cellType),'fontsize',14);
            
        end
    case{'movielinear','linearmovie'}

        resp = obj.mosaic{mosaicIndex}.get('response linear');
        ieMovie(resp);
        
    case{'raster','rasterresponse'}
        % Plot spike rasterse
        
        timeStep = obj.timing*dt;
        % bindur = ir.mosaic{1}.get('dt');
        
        for cellTypeInd = mosaicIndices
            
            nCells = obj.mosaic{cellTypeInd}.get('mosaic size');
            nCellStart = [1 1]; nCellEnd = nCells;
            if length(cellIndices)>0 
                nCellStart = [cellIndices];
                nCellEnd = [cellIndices];
            end
            nTrials = obj.mosaic{1}.get('numbertrials');
            
            lastspiketime = obj.mosaic{1}.get('lastspiketime');
            cellCtr = 0;
            
            for xcell = nCellStart(1):nCellEnd(1)
                for ycell = nCellStart(2):nCellEnd(2)
                    cellCtr = cellCtr+1;
                    if length(cellIndices)>2 
                        subplot(length(nCellStart(1):nCellEnd(1)),length(nCellStart(2):nCellEnd(2)),cellCtr);
                    end
                    for tr = 1:nTrials
                        % Get spike times
                        spikeTimes = obj.mosaic{cellTypeInd}.responseSpikes{xcell,ycell,tr};
                        hold on; 
                        
                        % Draw raster plots
                        % line([spikeTimes, spikeTimes]*timeStep,[tr tr-1],'color','k');                        
                        scatter([spikeTimes].*timeStep,[tr*ones(length(spikeTimes),1)],8,'o',color,'filled');
                        
                        axis([0 timeStep*lastspiketime 0 nTrials+1]);
                        
                    end
                    
                    if cellCtr == 1
                        
                        xlabel(sprintf('Time (sec)'));
                        ylabel(sprintf('Trial'));
                        title(sprintf('%s cell [%d %d] raster',obj.mosaic{cellTypeInd}.cellType,xcell,ycell));
                    end
                    
                end
            end
        end
        
    case{'psth'}
        timeStep = obj.timing;
                
        for cellTypeInd = mosaicIndices

        nCells = obj.mosaic{cellTypeInd}.get('mosaic size');
            nCellStart = [1 1]; nCellEnd = nCells;
            if length(cellIndices)>0 
                nCellStart = [cellIndices];
                nCellEnd = [cellIndices];
            end
            nTrials = obj.mosaic{1}.get('numbertrials');
            
            lastspiketime = obj.mosaic{1}.get('lastspiketime');
            cellCtr = 0;
            
            for xcell = nCellStart(1):nCellEnd(1)
                for ycell = nCellStart(2):nCellEnd(2)
                    cellCtr = cellCtr+1;
                    
                    psth = obj.mosaic{cellTypeInd}.get('psth');
                    
                    resp = RGB2XWFormat(psth(xcell,ycell,:));
                    
                    plot(timeStep:timeStep:timeStep*size(resp,2),resp);
                    
                    
                    if cellCtr == 1
                        xlabel(sprintf('Time (sec)'));
                        ylabel(sprintf('Conditional Intensity'));
                        title(sprintf('%s',obj.mosaic{cellTypeInd}.cellType));
                    end
                end
            end
        end
end
        
        