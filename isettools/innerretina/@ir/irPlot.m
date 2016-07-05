function irPlot(obj, varargin)
% Plots properties of the inner retina object
%
% Inputs: ir object, name of property to be plotted
%
% Outputs: plot(s)
%
% Properties that can be plotted:
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
%
%   osI = osCreate('identity');
%   innerRetina = irCreate(osI);
%   innerRetina.mosaicCreate('model','glm','type', 'on parasol');
%   innerRetina.compute(osI);
%
%   irPlot(innerRetina,'mosaic');
%   irPlot(innerRetina,'psth');
%   irPlot(innerRetina,'psth','type','onParasol');
%   irPlot(innerRetina,'psth','cell',[1 1]);
%   irPlot(innerRetina,'psth','type','onParasol','cell',[1 1]);
%
% (c) isetbio
% 09/2015 JRG
%% Parse input
p = inputParser;
p.CaseSensitive = false;

% Make key properties that can be set required arguments, and require
% values along with key names.
allowableFields = {...
    'rf',...
    'rfImage',...
    'mosaic',...
    'sRFcenter',...
    'sRFsurround',...
    'temporal','temporalfilter',...
    'ecc',...
    'tCenter',...
    'tSurround',...
    'postSpikeFilter',...
    'couplingFilter',...
    'responselinear','linear',...
    'nlResponse','nl','nonlinear','responsenl',...
    'responseSpikes','spike','voltage','responseVoltage',...
    'rasterResponse','raster','responseraster',...
    'psthResponse','psth','responsepsth'...
    };

% Special case ... please describe 'what' ...
p.addRequired('what',@(x) any(validatestring(ieParamFormat(x),allowableFields)));

p.addParameter('type',[],@ischar);      % Mosaic type?
p.addParameter('cell',[],@isnumeric);   % Vector indicating which cell
p.addParameter('hold','off',@ischar);   % Plotting related
p.addParameter('color','r',@ischar);

p.parse(varargin{:}); params = p.Results;
color      = p.Results.color; % params.color;
mosaicType = p.Results.type;  % params.type;
holdVal    = p.Results.hold;

% Parse mosaic type to allow for entering a number or a string
if ~isempty(mosaicType) && ischar(mosaicType)
    for cellTypeInd = 1:length(obj.mosaic)
        if strcmp(ieParamFormat(obj.mosaic{cellTypeInd}.cellType),ieParamFormat(mosaicType))
            mosaicTypeInd = cellTypeInd;
        end
    end
elseif ~isempty(mosaicType) && isnumeric(mosaicType)
    mosaicTypeInd = mosaicType;
else
    mosaicTypeInd = [];
end

cellID = params.cell;


if strcmpi(holdVal,'off')
    vcNewGraphWin([],'upperleftbig');
else
    hold on;
end


%% Set key-value pairs
switch ieParamFormat(params.what)
    case{'ecc'}
        %%% The temporal equivalent eccentricity of the retinal patch
        plotPatchEccentricity(obj.eyeAngle, obj.eyeRadius, obj.eyeSide, obj.temporalEquivEcc)
        
    case{'mosaic'}
        %%% Plot the mosaic of each RGC type
        
        % vcNewGraphWin([],'upperleftbig');
        % set(gcf,'position',[1000  540 893  798]);
        cmap = parula(16);
        
        for cellTypeInd = 1:length(obj.mosaic)
            
            spatialRFcontours = plotContours(obj.mosaic{cellTypeInd}, obj.spacing, obj.col);
            
            if length(obj.mosaic)>1
                subplot(ceil(length(obj.mosaic)/2),2,cellTypeInd);
            end
            
            nCells = size(obj.mosaic{cellTypeInd}.cellLocation);
            
            
            patchSizeX = obj.spacing;
            sensorRows = obj.col;
            
            umPerSensorPx = patchSizeX/sensorRows;
            
            for xcell = 1:nCells(1)
                for ycell = 1:nCells(2)
                    hold on;
                    % center
                    plot(umPerSensorPx*spatialRFcontours{xcell,ycell,1}(1,2:end),...
                        umPerSensorPx*spatialRFcontours{xcell,ycell,1}(2,2:end),...
                        'color',cmap(cellTypeInd,:));
                    %                     plot(umPerSensorPx*spatialRFcontours{xcell,ycell,1}(1,2:end),...
                    %                         umPerSensorPx*spatialRFcontours{xcell,ycell,1}(2,2:end),...
                    %                         'color','r','linewidth',5);
                    
                    %                     fill(umPerSensorPx*spatialRFcontours{xcell,ycell,1}(1,2:end),...
                    %                         umPerSensorPx*spatialRFcontours{xcell,ycell,1}(2,2:end),cmap(cellTypeInd,:));
                    hold on;
                    % surround
                    plot(umPerSensorPx*spatialRFcontours{xcell,ycell,2}(1,2:end),...
                        umPerSensorPx*spatialRFcontours{xcell,ycell,2}(2,2:end),...
                        'color',cmap(cellTypeInd+8,:));
                end
            end
            axis equal
            title(sprintf('%s',obj.mosaic{cellTypeInd}.cellType),'fontsize',14);
            xlabel(sprintf('Distance (\\mum)'),'fontsize',14);
            ylabel(sprintf('Distance (\\mum)'),'fontsize',14);
            xlabel('Distance (m)'); ylabel('Distance (m)'); set(gca,'fontsize',14);
            %             hl=legend('center','surround','location','ne'); set(hl,'fontsize',14);
        end
        % plot the cone mosaic with RGC contours!
        %         cone_mosaic = absorptions.human.coneType;
        %         [xg yg] = meshgrid([1:123,1:150]);
        %         xg2 = xg(1:123,1:150); yg2 = yg(1:123,1:150);
        %
        %         figure; scatter(xg2(:),yg2(:),40,4-cone_mosaic(:),'o','filled'); colormap jet; set(gca,'color',[0 0 0])
        %
        %         [xg yg] = meshgrid([1:90]); figure; scatter(xg(:),yg(:),40,4-cone_mosaic(:),'o','filled'); colormap jet; set(gca,'color',[0 0 0])
    case{'rf'}
        %%% A surface representing the RF (center - surround)
        % vcNewGraphWin([],'upperleftbig');
        % set(gcf,'position',[1000  540 893  798]);
        for cellTypeInd = 1%:length(obj.mosaic)
            
            %             subplot(3,2,cellTypeInd);
            surface(obj.mosaic{cellTypeInd}.sRFcenter{1,1}-obj.mosaic{cellTypeInd}.sRFsurround{1,1}); shading flat; view(40,40);
            title(sprintf('Spatial Receptive Field, %s',obj.mosaic{cellTypeInd}.cellType),'fontsize',14);
            xlabel(sprintf('Distance (\\mum)'),'fontsize',14);
            ylabel(sprintf('Distance (\\mum)'),'fontsize',14);
            zlabel(sprintf('Response (spikes/sec)'),'fontsize',14);
            axis([0 size(obj.mosaic{cellTypeInd}.sRFcenter{1,1},1) 0 size(obj.mosaic{cellTypeInd}.sRFcenter{1,1},2) -max(obj.mosaic{cellTypeInd}.sRFsurround{1,1}(:)) max(obj.mosaic{cellTypeInd}.sRFcenter{1,1}(:)) ]);
        end
        
    case{'rfimage'}
        % An image representing the RF surround
        % vcNewGraphWin([],'upperleftbig');
        % set(gcf,'position',[1000  540 893  798]);
        if ~isempty(mosaicTypeInd)
            cellTypeStart = mosaicTypeInd;
            cellTypeEnd = mosaicTypeInd;
        else
            cellTypeStart = 1;
            cellTypeEnd = length(obj.mosaic);
        end
        for cellTypeInd = cellTypeStart:cellTypeEnd
            
            if length(cellTypeStart:cellTypeEnd)>1
                subplot(ceil(length(cellTypeStart:cellTypeEnd)/2),2,cellTypeInd);
            end
            imagesc(obj.mosaic{cellTypeInd}.sRFcenter{1,1}-obj.mosaic{cellTypeInd}.sRFsurround{1,1});
            title(sprintf('Spatial Receptive Field, %s',obj.mosaic{cellTypeInd}.cellType),'fontsize',14);
            xlabel(sprintf('Distance (\\mum)'),'fontsize',14);
            ylabel(sprintf('Distance (\\mum)'),'fontsize',14);
            h = colorbar; ylabel(h, 'Response (spikes/sec)','fontsize',14);
            axis([0 size(obj.mosaic{cellTypeInd}.sRFcenter{1,1},1) 0 size(obj.mosaic{cellTypeInd}.sRFcenter{1,1},2)]);
        end
        
        
    case{'srfcenter'}
        
        %%% A surface representing the RF center
        % vcNewGraphWin([],'upperleftbig');
        % set(gcf,'position',[1000  540 893  798]);
        if ~isempty(mosaicTypeInd)
            cellTypeStart = mosaicTypeInd;
            cellTypeEnd   = mosaicTypeInd;
        else
            cellTypeStart = 1;
            cellTypeEnd = length(obj.mosaic);
        end
        
        for cellTypeInd = cellTypeStart:cellTypeEnd
            if length(cellTypeStart:cellTypeEnd)>1
                subplot(ceil(length(cellTypeStart:cellTypeEnd)/2),2,cellTypeInd);
            end
            
            % Let's fix this cellID logic.
            if ~isempty(cellID)
                nCells = cellID;
                xcellstart = cellID(1);
                ycellstart = cellID(2);
            else
                nCells = size(obj.mosaic{cellTypeInd}.cellLocation);
                xcellstart = 1; ycellstart = 1;
            end
            
            for xcell = xcellstart:nCells(1)
                for ycell = ycellstart:nCells(2)
                    
                    surface(obj.mosaic{cellTypeInd}.sRFcenter{xcell,ycell}); shading flat; view(40,40);
                    title(sprintf('Spatial Receptive Field, %s',obj.mosaic{cellTypeInd}.cellType),'fontsize',14);
                    xlabel(sprintf('Distance (\\mum)'),'fontsize',14);
                    ylabel(sprintf('Distance (\\mum)'),'fontsize',14);
                    zlabel(sprintf('Response (spikes/sec)'),'fontsize',14);
                    axis([0 size(obj.mosaic{cellTypeInd}.sRFcenter{1,1},1) 0 size(obj.mosaic{cellTypeInd}.sRFcenter{1,1},2) min(obj.mosaic{cellTypeInd}.sRFsurround{1,1}(:)) max(obj.mosaic{cellTypeInd}.sRFcenter{1,1}(:)) ]);
                    
                end
            end
        end
        
    case{'srfsurround'}
        %%% A surface representing the RF surround
        % vcNewGraphWin([],'upperleftbig');
        % set(gcf,'position',[1000  540 893  798]);
        if ~isempty(mosaicTypeInd)
            cellTypeStart = mosaicTypeInd;
            cellTypeEnd = mosaicTypeInd;
        else
            cellTypeStart = 1;
            cellTypeEnd = length(obj.mosaic);
        end
        for cellTypeInd = cellTypeStart:cellTypeEnd
            
            
            if length(cellTypeStart:cellTypeEnd)>1
                subplot(ceil(length(cellTypeStart:cellTypeEnd)/2),2,cellTypeInd);
            end
            surface(-obj.mosaic{cellTypeInd}.sRFsurround{1,1}); shading flat; view(40,40);
            title(sprintf('Spatial Receptive Field, %s',obj.mosaic{cellTypeInd}.cellType),'fontsize',14);
            xlabel(sprintf('Distance (\\mum)'),'fontsize',14);
            ylabel(sprintf('Distance (\\mum)'),'fontsize',14);
            zlabel(sprintf('Response (spikes/sec)'),'fontsize',14);
            axis([0 size(obj.mosaic{cellTypeInd}.sRFsurround{1,1},1) 0 size(obj.mosaic{cellTypeInd}.sRFsurround{1,1},2) min(-obj.mosaic{cellTypeInd}.sRFsurround{1,1}(:)) max(-obj.mosaic{cellTypeInd}.sRFsurround{1,1}(:)) ]);
        end
    case{'temporal','temporalfilter'}
        %%% Plot the RGB impulse response of each mosaic
        % vcNewGraphWin([],'upperleftbig');
        % set(gcf,'position',[1000  540 893  798]);
        if ~isempty(mosaicTypeInd)
            cellTypeStart = mosaicTypeInd;
            cellTypeEnd = mosaicTypeInd;
        else
            cellTypeStart = 1;
            cellTypeEnd = length(obj.mosaic);
        end
        for cellTypeInd = cellTypeStart:cellTypeEnd
            
            if length(cellTypeStart:cellTypeEnd)>1
                subplot(ceil(length(cellTypeStart:cellTypeEnd)/2),2,cellTypeInd);
            end
            %             plot((.01:.01:.2)-.2,bsxfun(@plus,horzcat(obj.mosaic{cellTypeInd}.tCenter{:}),[0 0 0.01]))
            cind = 'bgr';
            offset = [0 0 .01];
            hold on;
            for rgbInd = 1:3
                plot((.01:.01:.2)-.2,((obj.mosaic{cellTypeInd}.tCenter{rgbInd})+offset(rgbInd)),cind(rgbInd))
            end
            title(sprintf('Temporal Filter, RGB, %s',obj.mosaic{cellTypeInd}.cellType),'fontsize',14);
            xlabel(sprintf('Time (sec)'),'fontsize',14);
            ylabel(sprintf('Sensitivity (AU)'),'fontsize',14);
            legend('B','G','R','location','northwest');
        end
        
    case{'tcenter'}
        
        %%% Plot the RGB impulse response of each mosaic
        % vcNewGraphWin([],'upperleftbig');
        % set(gcf,'position',[1000  540 893  798]);
        
        timeStep = obj.timing;
        
        if ~isempty(mosaicTypeInd)
            cellTypeStart = mosaicTypeInd;
            cellTypeEnd = mosaicTypeInd;
        else
            cellTypeStart = 1;
            cellTypeEnd = length(obj.mosaic);
        end
        for cellTypeInd = cellTypeStart:cellTypeEnd
            
            if length(cellTypeStart:cellTypeEnd)>1
                subplot(ceil(length(cellTypeStart:cellTypeEnd)/2),2,cellTypeInd);
            end
            
            if ~isempty(cellID)
                nCells = cellID;
                xcellstart = cellID(1); ycellstart = cellID(2);
            else
                nCells = size(obj.mosaic{cellTypeInd}.cellLocation);
                xcellstart = 1; ycellstart = 1;
            end
            
            if strcmpi(ieParamFormat(obj.mosaic{cellTypeInd}.cellType),'offparasol')
                offMult = -1;
            else
                offMult = 1;
            end
            
            for xcell = xcellstart:nCells(1)
                for ycell = ycellstart:nCells(2)
                    if strcmpi(class(obj.mosaic{cellTypeInd}),'rgcphys')
                        plot(0:timeStep:timeStep*(-1+length(obj.mosaic{cellTypeInd}.tCenter{xcell,ycell})),offMult*((obj.mosaic{cellTypeInd}.tCenter{xcell,ycell})),'r','linewidth',4)
                        line([0 timeStep*(-1+length(obj.mosaic{cellTypeInd}.tCenter{xcell,ycell}))], [0 0]);
                    else
                        plot(0:timeStep:timeStep*(-1+length(obj.mosaic{cellTypeInd}.tCenter{1})),offMult*bsxfun(@plus,horzcat(obj.mosaic{cellTypeInd}.tCenter{:}),[0 0 0.01]))
                        line([0 timeStep*(-1+length(obj.mosaic{cellTypeInd}.tCenter{1}))], [0 0]);
                    end
                    title(sprintf('Temporal Impulse Response, RGB, %s',obj.mosaic{cellTypeInd}.cellType),'fontsize',14);
                    xlabel(sprintf('Time (sec)'),'fontsize',14);
                    ylabel(sprintf('Sensitivity (AU)'),'fontsize',14);
                end
            end
        end
        
    case{'tsurround'}
        
        timeStep = obj.timing;
        %%% Plot the RGB impulse response of each mosaic
        % vcNewGraphWin([],'upperleftbig');
        % set(gcf,'position',[1000  540 893  798]);
        if ~isempty(mosaicTypeInd)
            cellTypeStart = mosaicTypeInd;
            cellTypeEnd = mosaicTypeInd;
        else
            cellTypeStart = 1;
            cellTypeEnd = length(obj.mosaic);
        end
        for cellTypeInd = cellTypeStart:cellTypeEnd
            if length(cellTypeStart:cellTypeEnd)>1
                subplot(ceil(length(cellTypeStart:cellTypeEnd)/2),2,cellTypeInd);
            end
            
            if ~isempty(cellID)
                nCells = cellID;
                xcellstart = cellID(1); ycellstart = cellID(2);
            else
                nCells = size(obj.mosaic{cellTypeInd}.cellLocation);
                xcellstart = 1; ycellstart = 1;
            end
            
            for xcell = xcellstart:nCells(1)
                for ycell = ycellstart:nCells(2)
                    if strcmpi(class(obj.mosaic{cellTypeInd}),'rgcphys')
                        plot(0:timeStep:timeStep*(-1+length(obj.mosaic{cellTypeInd}.tSurround{xcell,ycell})),((obj.mosaic{cellTypeInd}.tSurround{xcell,ycell})),'r','linewidth',4)
                        line([0 timeStep*(-1+length(obj.mosaic{cellTypeInd}.tCenter{xcell,ycell}))], [0 0]);
                    else
                        plot(0:timeStep:timeStep*(-1+length(obj.mosaic{cellTypeInd}.tSurround{1})),bsxfun(@plus,horzcat(obj.mosaic{cellTypeInd}.tSurround{:}),[0 0 0.01]))
                        line([0 timeStep*(-1+length(obj.mosaic{cellTypeInd}.tCenter{1}))], [0 0]);
                    end
                    title(sprintf('Temporal Impulse Response, RGB, %s',obj.mosaic{cellTypeInd}.cellType),'fontsize',14);
                    xlabel(sprintf('Time (sec)'),'fontsize',14);
                    ylabel(sprintf('Sensitivity (AU)'),'fontsize',14);
                end
            end
        end
        
        
    case{'postspikefilter'}
        
        timeStep = .1*obj.timing;
        
        %%% Plot the post spike filter of each cell type
        % vcNewGraphWin([],'upperleftbig');
        % set(gcf,'position',[1000  540 893  798]);
        if ~isempty(mosaicTypeInd)
            cellTypeStart = mosaicTypeInd;
            cellTypeEnd = mosaicTypeInd;
        else
            cellTypeStart = 1;
            cellTypeEnd = length(obj.mosaic);
        end
        for cellTypeInd = cellTypeStart:cellTypeEnd
            
            if ~isempty(cellID)
                nCells = cellID;
                xcellstart = cellID(1); ycellstart = cellID(2);
            else
                nCells = size(obj.mosaic{cellTypeInd}.cellLocation);
                xcellstart = 1; ycellstart = 1;
            end
            
            for xcell = xcellstart:nCells(1)
                for ycell = ycellstart:nCells(2)
                    if  strcmpi(class(obj.mosaic{cellTypeInd}),'rgcphys')
                        psf =exp( squeeze(obj.mosaic{cellTypeInd}.postSpikeFilter{xcell,ycell}));
                    else
                        psf = squeeze(obj.mosaic{cellTypeInd}.couplingFilter{1,1}(1,1,:));
                    end
                    
                    if length(cellTypeStart:cellTypeEnd)>1
                        subplot(ceil(length(cellTypeStart:cellTypeEnd)/2),2,cellTypeInd);
                    end
                    plot(0:timeStep:timeStep*((-1+length(psf))), psf,'r','linewidth',4);
                    line([0 (-1+length(psf)).*timeStep],[ 1 1]);
                    title(sprintf('Exponentiated Post-Spike Filter, %s',obj.mosaic{cellTypeInd}.cellType),'fontsize',14);
                    xlabel(sprintf('Time (sec)'),'fontsize',14);
                    ylabel(sprintf('Gain (AU)'),'fontsize',14);
                end
            end
        end
        
    case{'couplingfilter'}
        
        
        %%% Plot the post spike filter of each cell type
        % vcNewGraphWin([],'upperleftbig');
        % set(gcf,'position',[1000  540 893  798]);
        if ~isempty(mosaicTypeInd)
            cellTypeStart = mosaicTypeInd;
            cellTypeEnd = mosaicTypeInd;
        else
            cellTypeStart = 1;
            cellTypeEnd = length(obj.mosaic);
        end
        for cellTypeInd = cellTypeStart:cellTypeEnd
            
            cplf = (horzcat(obj.mosaic{cellTypeInd}.couplingFilter{:}));
            
            if length(cellTypeStart:cellTypeEnd)>1
                subplot(ceil(length(cellTypeStart:cellTypeEnd)/2),2,cellTypeInd);
            end
            plot((1:length(cplf))./1000, squeeze(cplf(1,:,:)));
            title(sprintf('Exponentiated Coupling Filters, %s',obj.mosaic{cellTypeInd}.cellType),'fontsize',14);
            xlabel(sprintf('Time (sec)'),'fontsize',14);
            ylabel(sprintf('Response (spikes/sec)'),'fontsize',14);
            axis([0 0.6 -0.8 0.8]);
        end
        
    case{'responselinear','linear'}
        % vcNewGraphWin([],'upperleftbig');
        % set(gcf,'position',[1000  540 893  798]);
        
        timeStep = obj.timing;
        
        if ~isempty(mosaicTypeInd)
            cellTypeStart = mosaicTypeInd;
            cellTypeEnd = mosaicTypeInd;
        else
            cellTypeStart = 1;
            cellTypeEnd = length(obj.mosaic);
        end
        for cellTypeInd = cellTypeStart:cellTypeEnd
            if ~isempty(cellID)
                nCells = cellID;
                xcellstart = cellID(1); ycellstart = cellID(2);
            else
                nCells = size(obj.mosaic{cellTypeInd}.cellLocation);
                xcellstart = 1; ycellstart = 1;
            end
            for xcell = xcellstart:nCells(1)
                for ycell = ycellstart:nCells(2)
                    
                    meanVoltage{xcell,ycell} = ((obj.mosaic{cellTypeInd}.responseLinear{xcell,ycell}));
                end
            end
            if length(cellTypeStart:cellTypeEnd)>1
                subplot(ceil(length(cellTypeStart:cellTypeEnd)/2),2,cellTypeInd);
            end
            
            plot(timeStep:timeStep:timeStep*length(vertcat(meanVoltage{:})),vertcat(meanVoltage{:})');
            xlabel(sprintf('Time (msec)'),'fontsize',14);
            ylabel(sprintf('Membrane Voltage (\\muV)'),'fontsize',14);
            title(sprintf('%s',obj.mosaic{cellTypeInd}.cellType),'fontsize',14);
            
        end
        
        
    case{'nlresponse','nl','nonlinear','responsenl'}
        % vcNewGraphWin([],'upperleftbig');
        % set(gcf,'position',[1000  540 893  798]);
        
        if ~isempty(mosaicTypeInd)
            cellTypeStart = mosaicTypeInd;
            cellTypeEnd = mosaicTypeInd;
        else
            cellTypeStart = 1;
            cellTypeEnd = length(obj.mosaic);
        end
        for cellTypeInd = cellTypeStart:cellTypeEnd
            if ~isempty(cellID)
                nCells = cellID;
                xcellstart = cellID(1); ycellstart = cellID(2);
            else
                nCells = size(obj.mosaic{cellTypeInd}.cellLocation);
                xcellstart = 1; ycellstart = 1;
            end
            for xcell = xcellstart:nCells(1)
                for ycell = ycellstart:nCells(2)
                    
                    meanVoltage{xcell,ycell} = ((obj.mosaic{cellTypeInd}.nlResponse{xcell,ycell}));
                end
            end
            if length(cellTypeStart:cellTypeEnd)>1
                subplot(ceil(length(cellTypeStart:cellTypeEnd)/2),2,cellTypeInd);
            end
            plot(vertcat(meanVoltage{:})');
            % set(gca,'yscale','log');
            xlabel(sprintf('Time (msec)'),'fontsize',14);
            ylabel(sprintf('Membrane Voltage (\\muV)'),'fontsize',14);
            title(sprintf('%s',obj.mosaic{cellTypeInd}.cellType),'fontsize',14);
            
        end
        
    case{'responseSpikes','responsevoltage','spike','voltage'}
        
        %%% Plot the membrane voltages for a random trial
        
        % vcNewGraphWin([],'upperleftbig');
        % set(gcf,'position',[1000  540 893  798]);
        % szSpike = size(horzcat(obj.mosaic{1}.responseVoltage{1,1,:}));
        
        if ~isempty(mosaicTypeInd)
            cellTypeStart = mosaicTypeInd;
            cellTypeEnd = mosaicTypeInd;
        else
            cellTypeStart = 1;
            cellTypeEnd = length(obj.mosaic);
        end
        for cellTypeInd = cellTypeStart:cellTypeEnd
            clear meanVoltage
            if ~isempty(cellID)
                nCells = cellID;
                xcellstart = cellID(1); ycellstart = cellID(2);
            else
                nCells = size(obj.mosaic{cellTypeInd}.cellLocation);
                xcellstart = 1; ycellstart = 1;
            end
            
            meanVoltage = cell(length(xcellstart:nCells(1)),length(ycellstart:nCells(2)));
            for xcell = xcellstart:nCells(1)
                for ycell = ycellstart:nCells(2)
                    meanVoltage{xcell,ycell} = ((obj.mosaic{cellTypeInd}.responseVoltage{xcell,ycell}));
                end
            end
            
            if length(cellTypeStart:cellTypeEnd)>1
                subplot(ceil(length(cellTypeStart:cellTypeEnd)/2),2,cellTypeInd);
            end
            mv = horzcat(meanVoltage{:});
            plot((1:length(mv))/100,mv);
            
            xlabel(sprintf('Time (msec)'),'fontsize',14);
            ylabel(sprintf('exp(Membrane Voltage (\\muV))'),'fontsize',14);
            title(sprintf('%s',obj.mosaic{cellTypeInd}.cellType));%,'fontsize',14);
            
            clear meanVoltage
        end
        
    case{'raster'}
        % There may be some dt issue to deal with here.
        % @JRG. This section is way too long and needs to be clarified and
        % functionalized.
        
        if ~isempty(mosaicTypeInd)
            cellTypeStart = mosaicTypeInd;
            cellTypeEnd = mosaicTypeInd;
        else
            cellTypeStart = 1;
            cellTypeEnd = length(obj.mosaic);
        end
        
        % Let's just do one mosaic.
        for cellTypeInd = cellTypeStart:cellTypeEnd
            
            if strcmpi(class(obj.mosaic{cellTypeInd}),'rgcphys');
                bindur = .01/1.208;
            else
                bindur = .01;
            end
            
            numberTrials = obj.mosaic{cellTypeInd}.get('numberTrials');
            
            % Why this?
            set(gcf,'position',[0.0069    0.3933    0.8021    0.5000]);
            cellCtr = 0;
            
            clear psth tsp mtsp
            clear meanVoltage
            maxTrials = size(obj.mosaic{cellTypeInd}.responseSpikes,3);
            
            if ~isempty(cellID)
                nCells = cellID;
                xcellstart = cellID(1); ycellstart = cellID(2);
            else
                nCells = size(obj.mosaic{cellTypeInd}.responseSpikes);
                % Sometimes @JRG adjusts this.  Some EJ data thing, I
                % suppose.  Maybe we can get rid of it here.
                if ~strcmpi(class(obj.mosaic{cellTypeInd}),'rgcphys') && ...
                        nCells(2) == 1;
                    nCells(1) = ceil(sqrt(nCells(1)));
                    nCells(2) = nCells(1);
                end
                xcellstart = 1; ycellstart = 1;
            end
            
            for xcell = xcellstart:nCells(1)
                for ycell = ycellstart:nCells(2)
                    cellCtr = cellCtr+1;
                    [jv,iv] = ind2sub([nCells(1),nCells(2)],cellCtr);
                    cellCtr2 = sub2ind([nCells(2),nCells(1)],iv,jv);
                    
                    % For each trial
                    for tr = 1:numberTrials;
                        clear spikeTimesP
                        
                        if strcmpi(holdVal,'off')
                            subplot(length(ycellstart:nCells(2)),length(xcellstart:nCells(1)),cellCtr2);
                        end
                        
                        if strcmpi(class(obj.mosaic{cellTypeInd}),'rgcphys');
                            % @JRG:  Why is there a 0.1 here?
                            spikeTimesP = .1*(obj.mosaic{cellTypeInd}.responseSpikes{xcell,ycell,tr,1});
                        else
                            spikeTimesP = (obj.mosaic{cellTypeInd}.responseSpikes{xcell,ycell,tr,1});
                        end
                        
                        % @JRG.  ???
                        if length(spikeTimesP) == 2
                            spikeTimesP = [spikeTimesP; 0];
                        end
                        
                        if ~isempty(spikeTimesP)
                            hold on;
                            scatter(spikeTimesP.*bindur,tr*ones(length(spikeTimesP),1),8,'o',color,'filled');
                        end
                        
                        if xcell == 1 && ycell == 1
                            xlabel('Time (sec)'); ylabel('Trial');
                            title(sprintf('%s cell [%d %d]',obj.mosaic{cellTypeInd}.cellType,xcell,ycell));
                        end
                    end
                    
                    if isempty(obj.mosaic{cellTypeInd}.responseSpikes)
                        maxt = length((obj.mosaic{cellTypeInd}.responseVoltage{1,1}));
                        axis([0 .00001*maxt 0 maxTrials]);
                    else
                        maxt = max((obj.mosaic{cellTypeInd}.responseSpikes{1}));
                        axis([0 .01*maxt 0 maxTrials]);
                    end
                end
            end
        end
        
        
    case{'psth'}
        % Post-stimulus time histogram
        % @JRG.  Broken now.
        %
        % Example:
        %   irPlot(ir,'psth','type','on midget');
        
        % @JRG:  I am forcing this to 1 for now.  We need to be able to
        % pass it in.
        cellTypeInd = 1;
        
        nCells = obj.mosaic{cellTypeInd}.get('mosaic size');
        % Which of the cells.  I think there should only be one class.
        psth = obj.mosaic{cellTypeInd}.get('psth');
        thisPlot = 1;
        t = size(psth,3);
        hf = vcNewGraphWin;
        for ii=1:nCells(1)
            for jj=1:nCells(2)
                subplot(nCells(1), nCells(2),thisPlot)
                spikes = squeeze(psth(ii,jj,:));
                plot(1:t,spikes)
                if thisPlot == 1
                    xlabel('Time (sec)','FontSize',14); ylabel(sprintf('PSTH\n(spikes/sec)'),'FontSize',14);
                    set(hf,'Name',sprintf('%s cell type',obj.mosaic{cellTypeInd}.cellType));
                end
                thisPlot = thisPlot + 1;
                set(gca,'xlim',[0 t],'ylim',[0 max(spikes)],'ytick',[0 max(spikes)],'xtick',[0 t]);
                set(gca,'fontsize',12)
            end
        end
end

end


