function [uData, hf] = plot(obj, plotType, varargin)
% Plot function for rgcMosaic
%    [uData, hf] = rgcMosaic.plot(plotType, varargin)
%
% Required input
%   plotType - string, type of plot

% Optional input (key-val pairs in varargin):
%   'hf' - figure or axis handle
%      By default, this is the figure handle of the rgcMosaic window
%      It can also be one of the two axes in the window.
%      Otherwise
%             []            - create plot in new figure
%             figure handle - an alternative figure
%             'none'        - don't plot, only generate uData
%   'gam' - Display gamma parameter (imagesc(x.^gam));
%
% Outputs:
%   uData - Computed user data
%   hf    - figure handle
%
% Plot type can be chosen from
%   'spike mean image'    - Gray scale image of mean spikes per
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
p.addParameter('hf', obj.fig, @isgraphics);  % figure handle
p.addParameter('gam',1,@isnumeric);

p.parse(plotType, varargin{:});
hf  = p.Results.hf;
gam = p.Results.gam;

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
        spikes = obj.get('spikes');  % Number of spikes in each ms
        if isempty(spikes)
            disp('No spikes have been computed (responseSpikes missing)');
            return;
        end
        
        img = mean(spikes,3);  % Mean spikes per millisecond
        img = img*1000;        % Mean spikes per second
        colormap(gray(256)); imagesc(img.^gam); axis image;
        axis off; colorbar; drawnow;
        title('Spikes/sec');
        
    case 'spikemovie'
        % Movie of spiking activity
        % Should this work with the mosaic plot of circles?
        psthTest = obj.get('spikes');
       
        clear vParams; vParams = [];
        vParams.gamma     = gam;
        vParams.FrameRate = 30; 
        vParams.show      = true; %vParams.step = 2;
        frameSkip = round(1./obj.get('dt'));
        ieMovie(psthTest(:,:,1:frameSkip:end),vParams);
        
    case{'linearmovie'}
        % Continuous voltages prior to spike generation
        responseLinear = obj.get('responseLinear');
        if isempty(responseLinear)
            disp('No response computed');
            return;
        end
        
        clear vParams;
        vParams.FrameRate = 30;
        vParams.show      = true;
        vParams.gamma     = gam;

        % Should we add controls like the cone mosaic window, or just play
        % it once?
        ieMovie(responseLinear,vParams);
        
    case 'psthmeanimage'
        psth = obj.get('psth');
        imagesc(mean(psth,3).^gam);
        axis image; colormap(gray(256)); colorbar;
        set(gca,'xticklabels','','yticklabels','');
        drawnow;
        
    case 'psth'
        % Peri-stimulus time graph for all cells.
        % Kind of a weird plot to make.
        % Need a flag to force a new window;
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
        % Plot the mosaic spatial receptive field geometry, with no fill.
        % We need to add the possibility of elliptical forms some day.
        % And we should figure out how to do center/surround
        cla reset;
        
        % We should subsample the number shown when there are many.
        center = obj.cellLocation;  % um w.r.t. center of image
        
        % Why are we plotting sqrt(2)* radius?  1 std is radius.
        % But 1.4 radius looks good.  Is that why?
        % (BW)
        radius = obj.rfDiameter/2;
        ellipseMatrix = obj.ellipseMatrix;
        ieShape('ellipse','center',center,...
            'radius',sqrt(2)*radius,...
            'ellipseParameters',vertcat(ellipseMatrix{:}),...
            'color','b');
        
        % Sets the axis limits
        set(gca,...
            'xlim',[min(center(:,2)) - 3*radius, max(center(:,2)) + 3*radius],...
            'ylim',[min(center(:,1)) - 3*radius, max(center(:,1)) + 3*radius]);
        xlabel(sprintf('Distance (\\mum)'),'fontsize',14);
        
    case 'mosaicfill'
        % Plot the mosaic spatial receptive field geometry but shading the
        % interior to reflect the activity level.
        % Perhaps we should take in a 'fill' parameter, rather than have
        % two different case statements.
        cla reset;
        
        spikes = obj.get('spikes');  % Number of spikes in each ms
        if isempty(spikes)
            disp('No spikes have been computed.');
            return;
        end
        img = mean(spikes,3);  % Mean spikes per millisecond
        img = img*1000;        % Mean spikes per second
        img = img.^gam;        % This makes the units arbitrary
        
        umPerSample = 1e6*obj.get('meters per sample');
        
        % We should subsample the number shown when there are many.
        center = RGB2XWFormat(obj.cellLocation);  % um w.r.t. center of image
        center = center * diag(umPerSample);
        
        radius = obj.rfDiameter/2;
        radius = radius * mean(umPerSample);
        
        ellipseMatrix = obj.ellipseMatrix;
        ieShape('ellipse','center',center,...
            'radius',sqrt(2)*radius,...
            'ellipseParameters',vertcat(ellipseMatrix{:}),...
            'fillArray',flipud(img)); %obj.responseLinear(:,:,10)
        
        % Sets the plot axis limits
        set(gca,...
            'xlim',[min(center(:,2)) - 3*radius, max(center(:,2)) + 3*radius],...
            'ylim',[min(center(:,1)) - 3*radius, max(center(:,1)) + 3*radius]);
        xlabel(sprintf('Distance (\\mum)'),'fontsize',14);
        
        colormap(gray(256)); colorbar; drawnow;
        title('Mean spikes (a.u.)');
        
    case 'mosaicsurf'
        % Plots a sub-sampled set of spatial RF as surface (mesh) plots
        % Should make the skip parameters selectable
        
        % TO DEBUG with new cellLocation format.
        
        rfCoords = vertcat(obj.cellLocation);
        rfMinR = min(rfCoords(:,1)); rfMaxR = max(rfCoords(:,1));
        rfMinC = min(rfCoords(:,2)); rfMaxC = max(rfCoords(:,2));
        
        rfSize = size(obj.sRFcenter{1,1});
        
        edgePadding = 4;
        spStim = zeros(edgePadding+ceil(rfSize(1)/1)+ceil(rfMaxR-rfMinR),edgePadding+ceil(rfSize(2)/1)+ceil(rfMaxC-rfMinC));
        
        % Sub-sampling values
        startInd = 2; skipInd = 1;
        [row,col] = size(obj.cellLocation);
        
        % Initialize all the cell arrays to be the same
        rvStart = cell(row,col); rvEnd = rvStart;
        cvStart = rvStart; cvEnd = rvStart;
        
        for ri = startInd:skipInd:size(obj.cellLocation,1)
            for ci = startInd:skipInd:size(obj.cellLocation,2)
                %         [ri ci]
                rvStart{ri,ci} = 1+ceil(obj.cellLocation{ri,ci}(1) +ceil((rfMaxR-rfMinR)/2)+1);% - ceil(rfSize(1)/2)+1);
                rvEnd{ri,ci}   = 1+ceil(obj.cellLocation{ri,ci}(1) +ceil((rfMaxR-rfMinR)/2)) + ceil(rfSize(1)/1);
                
                cvStart{ri,ci} = 1+ceil(obj.cellLocation{ri,ci}(2) +ceil((rfMaxC-rfMinC)/2)+1);% - ceil(rfSize(2)/2)+1);
                cvEnd{ri,ci}   = 1+ceil(obj.cellLocation{ri,ci}(2) +ceil((rfMaxC-rfMinC)/2) + ceil(rfSize(2)/1));
                
                if (rvStart{ri,ci} > 0) && (cvStart{ri,ci} > 0)
                    
                    spStim(rvStart{ri,ci}:rvEnd{ri,ci},cvStart{ri,ci}:cvEnd{ri,ci}) = ...
                        spStim(rvStart{ri,ci}:rvEnd{ri,ci},cvStart{ri,ci}:cvEnd{ri,ci})+obj.sRFcenter{ri,ci}-obj.sRFsurround{ri,ci};
                end
            end
        end
        
        if skipInd ~= 1
            % Let the user know if you are skippiing plots
            fprintf('Skip index (skipInd) = %d\n',skipInd);
        end
        imagesc(spStim);
        
    otherwise
        error('Unknown plot type %s\n',plotType);
end

end


