function val = rgcPlot(obj, varargin)
% rgcMosaicGet: a method of @rgcMosaic that gets rgcMosaic object 
% parameters using the input parser structure.
% 
% Parameters:
%       {''} -  
% 
% 9/2015 JRG 

% Check for the number of arguments and create parser object.
% Parse key-value pairs.
% 
% Check key names with a case-insensitive string, errors in this code are
% attributed to this function and not the parser object.
error(nargchk(0, Inf, nargin));
p = inputParser; p.CaseSensitive = false; p.FunctionName = mfilename;

% This flag causes the parser not to throw an error here in the superclass
% call. The subclass call will throw an error.
% p.KeepUnmatched = true;

% Make key properties that can be set required arguments, and require
% values along with key names.
allowableFieldsToSet = {...
        'rf',...
        'mosaic',...
        'sRFcenter',...
        'sRFsurround',...
        'ir',...
        'tCenter',...
        'tSurround',...
        'postSpikeFilter',...
        'couplingFilter',...
        'linearResponse',...
        'nlResponse',...
        'spikeResponse',...
        'rasterResponse',...
        'psthResponse'...
    };
p.addRequired('what',@(x) any(validatestring(x,allowableFieldsToSet)));

% % Define what units are allowable.
% allowableUnitStrings = {'a', 'ma', 'ua', 'na', 'pa'}; % amps to picoamps
% 
% % Set up key value pairs.
% % Defaults units:
% p.addParameter('units','pa',@(x) any(validatestring(x,allowableUnitStrings)));

% Parse and put results into structure p.
p.parse(varargin{:}); params = p.Results;

% % Old error check on input.
% if ~exist('params','var') || isempty(params)
%     error('Parameter field required.');
% end
% if ~exist('val','var'),   error('Value field required.'); end;

% Set key-value pairs.
switch lower(params.what)
    
    case{'mosaic'}
        %%% Plot the mosaic of each RGC type
                
        figure;
        set(gcf,'position',[1000  540 893  798]);
        cmap = parula(16);
        
        for cellTypeInd = 1:length(obj.mosaic)
            
            spatialRFcontours = plotContours(obj.mosaic{cellTypeInd});
            
            subplot(3,2,cellTypeInd);
            
            nCells = size(obj.mosaic{cellTypeInd}.cellLocation);
            for xcell = 1:nCells(1)
                for ycell = 1:nCells(2)
                    hold on;
                    % center
                    plot(spatialRFcontours{xcell,ycell,1}(1,2:end),...
                        spatialRFcontours{xcell,ycell,1}(2,2:end),...
                        'color',cmap(cellTypeInd,:));
                    hold on;
                    % surround
                    plot(spatialRFcontours{xcell,ycell,2}(1,2:end),...
                        spatialRFcontours{xcell,ycell,2}(2,2:end),...
                        'color',cmap(cellTypeInd+8,:));
                end
            end
            axis equal
            title(sprintf('%s',obj.mosaic{cellTypeInd}.cellType),'fontsize',16);
            xlabel(sprintf('Distance (\\mum)'),'fontsize',16);
            ylabel(sprintf('Distance (\\mum)'),'fontsize',16);
        end
        
    case{'rf'}
        %%% A surface representing the RF surround
        figure;
        set(gcf,'position',[1000  540 893  798]);
        for cellTypeInd = 1:length(obj.mosaic)
            
            subplot(3,2,cellTypeInd);
            surface(obj.mosaic{cellTypeInd}.sRFcenter{1,1}-obj.mosaic{cellTypeInd}.sRFsurround{1,1}); shading flat; view(40,40);
            title(sprintf('Spatial Receptive Field, %s',obj.mosaic{cellTypeInd}.cellType),'fontsize',16);
            xlabel(sprintf('Distance (\\mum)'),'fontsize',16);
            ylabel(sprintf('Distance (\\mum)'),'fontsize',16);
            zlabel(sprintf('Response (spikes/sec)'),'fontsize',16);
            axis([0 size(obj.mosaic{5}.sRFcenter{1,1},1) 0 size(obj.mosaic{5}.sRFcenter{1,1},2) -max(obj.mosaic{cellTypeInd}.sRFsurround{1,1}(:)) max(obj.mosaic{cellTypeInd}.sRFcenter{1,1}(:)) ]);
        end
        
    case{'srfcenter'}
        
        %%% A surface representing the RF center
        figure;
        set(gcf,'position',[1000  540 893  798]);
        for cellTypeInd = 1:length(obj.mosaic)
            
            subplot(3,2,cellTypeInd); 
            surface(obj.mosaic{cellTypeInd}.sRFcenter{1,1}); shading flat; view(40,40);
            title(sprintf('Spatial Receptive Field, %s',obj.mosaic{cellTypeInd}.cellType),'fontsize',16);
            xlabel(sprintf('Distance (\\mum)'),'fontsize',16);
            ylabel(sprintf('Distance (\\mum)'),'fontsize',16);
            zlabel(sprintf('Response (spikes/sec)'),'fontsize',16);
            axis([0 size(obj.mosaic{5}.sRFcenter{1,1},1) 0 size(obj.mosaic{5}.sRFcenter{1,1},2) 0 max(obj.mosaic{cellTypeInd}.sRFcenter{1,1}(:)) ]);
        end

    case{'srfsurround'}
        %%% A surface representing the RF surround
        figure;
        set(gcf,'position',[1000  540 893  798]);
        for cellTypeInd = 1:length(obj.mosaic)
            
            subplot(3,2,cellTypeInd); 
            surface(obj.mosaic{cellTypeInd}.sRFsurround{1,1}); shading flat; view(40,40);
            title(sprintf('Spatial Receptive Field, %s',obj.mosaic{cellTypeInd}.cellType),'fontsize',16);
            xlabel(sprintf('Distance (\\mum)'),'fontsize',16);
            ylabel(sprintf('Distance (\\mum)'),'fontsize',16);
            zlabel(sprintf('Response (spikes/sec)'),'fontsize',16);
            axis([0 size(obj.mosaic{5}.sRFsurround{1,1},1) 0 size(obj.mosaic{5}.sRFsurround{1,1},2) 0 max(obj.mosaic{cellTypeInd}.sRFsurround{1,1}(:)) ]);
        end
    case{'ir'}        
        %%% Plot the RGB impulse response of each mosaic
        figure;     
        set(gcf,'position',[1000  540 893  798]);
        for cellTypeInd = 1:length(obj.mosaic)
            
            subplot(3,2,cellTypeInd);
            plot(.01:.01:.2,bsxfun(@plus,horzcat(obj.mosaic{cellTypeInd}.tCenter{:}),[0 0 0.01]))
            title(sprintf('Temporal Impulse Response, RGB, %s',obj.mosaic{cellTypeInd}.cellType),'fontsize',16);
            xlabel(sprintf('Time (sec)'),'fontsize',16);
            ylabel(sprintf('Response (spikes/sec)'),'fontsize',16);
        end
        
    case{'tcenter'}
        
        %%% Plot the RGB impulse response of each mosaic
        figure;
        set(gcf,'position',[1000  540 893  798]);        
        for cellTypeInd = 1:length(obj.mosaic)
            
            subplot(3,2,cellTypeInd);
            plot(.01:.01:.2,bsxfun(@plus,horzcat(obj.mosaic{cellTypeInd}.tCenter{:}),[0 0 0.01]))
            title(sprintf('Temporal Impulse Response, RGB, %s',obj.mosaic{cellTypeInd}.cellType),'fontsize',16);
            xlabel(sprintf('Time (sec)'),'fontsize',16);
            ylabel(sprintf('Response (spikes/sec)'),'fontsize',16);
        end
        
    case{'tsurround'}
        
        %%% Plot the RGB impulse response of each mosaic
        figure;   
        set(gcf,'position',[1000  540 893  798]);
        for cellTypeInd = 1:length(obj.mosaic)
            
            subplot(3,2,cellTypeInd);
            plot(.01:.01:.2,bsxfun(@plus,horzcat(obj.mosaic{cellTypeInd}.tSurround{:}),[0 0 0.01]))
            title(sprintf('Temporal Impulse Response, RGB, %s',obj.mosaic{cellTypeInd}.cellType),'fontsize',16);
            xlabel(sprintf('Time (sec)'),'fontsize',16);
            ylabel(sprintf('Response (spikes/sec)'),'fontsize',16);
        end
        
        
    case{'postspikefilter'}
        
        %%% Plot the post spike filter of each cell type
        figure;
        set(gcf,'position',[1000  540 893  798]);
        for cellTypeInd = 1:length(obj.mosaic)
            psf = squeeze(obj.mosaic{cellTypeInd}.couplingFilter{1,1}(1,1,:));
            
            subplot(3,2,cellTypeInd);
            plot((1:length(psf))./1000, psf);
            title(sprintf('Exponentiated Post-Spike Filter, %s',obj.mosaic{cellTypeInd}.cellType),'fontsize',16);
            xlabel(sprintf('Time (sec)'),'fontsize',16);
            ylabel(sprintf('Response (spikes/sec)'),'fontsize',16);
        end
        
    case{'couplingfilter'}
        
        
        %%% Plot the post spike filter of each cell type
        figure;
        set(gcf,'position',[1000  540 893  798]);
        for cellTypeInd = 1:length(obj.mosaic)
            
            cplf = (horzcat(obj.mosaic{cellTypeInd}.couplingFilter{:}));
            
            subplot(3,2,cellTypeInd);
            plot((1:length(cplf))./1000, squeeze(cplf(1,:,:)));
            title(sprintf('Exponentiated Coupling Filters, %s',obj.mosaic{cellTypeInd}.cellType),'fontsize',16);
            xlabel(sprintf('Time (sec)'),'fontsize',16);
            ylabel(sprintf('Response (spikes/sec)'),'fontsize',16);
            axis([0 0.6 -0.8 0.8]);
        end
        
    case{'linearresponse'}
        figure;
        set(gcf,'position',[1000  540 893  798]);
        
        for cellTypeInd = 1:length(obj.mosaic)
            nCells = size(obj.mosaic{cellTypeInd}.cellLocation);
            for xcell = 1:nCells(1)
                for ycell = 1:nCells(2)
                    
                    meanVoltage{xcell,ycell} = ((obj.mosaic{cellTypeInd}.linearResponse{xcell,ycell}));
                end
            end
            subplot(3,2,cellTypeInd);
            plot(vertcat(meanVoltage{:})');
            xlabel(sprintf('Time (msec)'),'fontsize',16);
            ylabel(sprintf('Membrane Voltage (\\muV)'),'fontsize',16);
            title(sprintf('%s',obj.mosaic{cellTypeInd}.cellType),'fontsize',16);
            
        end
        
                
    case{'nlresponse'}
        figure;
        set(gcf,'position',[1000  540 893  798]);
        
        for cellTypeInd = 1:length(obj.mosaic)
            nCells = size(obj.mosaic{cellTypeInd}.cellLocation);
            for xcell = 1:nCells(1)
                for ycell = 1:nCells(2)
                    
                    meanVoltage{xcell,ycell} = ((obj.mosaic{cellTypeInd}.nlResponse{xcell,ycell}));
                end
            end
            subplot(3,2,cellTypeInd);
            plot(vertcat(meanVoltage{:})');
            % set(gca,'yscale','log');
            xlabel(sprintf('Time (msec)'),'fontsize',16);
            ylabel(sprintf('Membrane Voltage (\\muV)'),'fontsize',16);
            title(sprintf('%s',obj.mosaic{cellTypeInd}.cellType),'fontsize',16);
            
        end
        
    case{'spikeresponse'}
        
        %%% Plot the membrane voltages for a random trial
        
        figure;        
        set(gcf,'position',[1000  540 893  798]);
        szSpike = size(horzcat(obj.mosaic{1}.spikeResponse{1,1,:,2}));
        
        for cellTypeInd = 1:length(obj.mosaic)
            clear meanVoltage
            nCells = size(obj.mosaic{cellTypeInd}.cellLocation);
            for xcell = 1:nCells(1)
                for ycell = 1:nCells(2)
                    
                    meanVoltage{xcell,ycell} = mean(horzcat(obj.mosaic{cellTypeInd}.spikeResponse{xcell,ycell,:,2}),2);
                end
            end
            subplot(3,2,cellTypeInd);
            plot(horzcat(meanVoltage{:}));
            xlabel(sprintf('Time (msec)'),'fontsize',16);
            ylabel(sprintf('Membrane Voltage (\\muV)'),'fontsize',16);
            title(sprintf('%s',obj.mosaic{cellTypeInd}.cellType),'fontsize',16);
            
        end
        
    case{'rasterresponse'}
        
        
        dt = .01; % make this a get from sensor
        
        for cellTypeInd = 1:length(obj.mosaic)
            
            figure;
            set(gcf,'position',[1000  540 893  798]);
            cellCtr = 0;
            clear psth tsp mtsp
            
            nCells = size(obj.mosaic{cellTypeInd}.cellLocation);
            szSpike = size(horzcat(obj.mosaic{1}.spikeResponse{1,1,:,2}));
            maxTrials = szSpike(2);
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
                    % [psth{xcell,ycell},tt,pstv,spr] = compPSTH(mtsp*dt, .001, .002, [0 1], .005);
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
            suptitle(sprintf('%s',obj.mosaic{cellTypeInd}.cellType));
        end
        
        
    case{'psthresponse'}
        
        
        dt = .01; % make this a get from sensor
        
        for cellTypeInd = 1:length(obj.mosaic)
            
            figure;
            set(gcf,'position',[1000  540 893  798]);
            cellCtr = 0;
            clear psth tsp mtsp
            
            nCells = size(obj.mosaic{cellTypeInd}.cellLocation);
            szSpike = size(horzcat(obj.mosaic{1}.spikeResponse{1,1,:,2}));
            maxTrials = szSpike(2);
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
                    % axis([0 30 0 maxTrials]);
                    
                    
                    % subplot(nCells(1),nCells(2),cellCtr);
                                plot(tt/.01,psth{xcell,ycell});
                                if ~isnan(psth{xcell,ycell})
                                    axis([0 30 0 max(psth{xcell,ycell})]);
                                end
                    
                    % hold on;
                    % gf = obj.mosaic{1}.generatorFunction;
                    % plot((obj.mosaic{cellTypeInd}.linearResponse{xcell,ycell}));
                    % axis([0 30 min(obj.mosaic{cellTypeInd}.linearResponse{xcell,ycell}) max(obj.mosaic{cellTypeInd}.linearResponse{xcell,ycell})]);
                    
                end
            end
            % mosaicSet(obj.mosaic{cellTypeInd},'rasterResponse',raster);
            % mosaicSet(obj.mosaic{cellTypeInd},'psthResponse',psth);
            suptitle(sprintf('%s',obj.mosaic{cellTypeInd}.cellType));
        end
end

return;

