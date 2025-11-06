function [uData, hdl] = plot(cmosaic,plotType, allE, varargin)
% cMosaic.plot - Plot methods for a cMosaic
%
%   [uData, hdl] = plot(cmosaic, plotType, allE, varargin)
%
% Description
%   Provides several visualization/analysis plots for a cMosaic object,
%   including activation maps, line profiles, ROI visualizations, and
%   spectral/photopigment functions. Unless otherwise noted, a figure
%   window is created (or reused when 'hdl' is supplied) and a struct of
%   useful data (uData) is returned.
%
% Inputs
%   cmosaic   (cMosaic)  A cMosaic instance.
%   plotType  (char)     Plot selector (see “Plot types” below).
%   allE      (numeric)  Precomputed excitations/activations array passed
%                        through to subroutines that require it. May be [] 
%                        for plot types that do not use it.
%
% Name-Value Pairs (varargin)
%   'roi'          : regionOfInterest object used by ROI plots/line cuts.
%   'conetype'     : char or cellstr of cone types to include: {'l','m','s'}
%                    (default {'l','m','s'}). Single types like 'l' allowed.
%   'hdl'          : matlab.ui.Figure handle to draw into (default [] → new).
%   'plottitle'    : char title string (default '').
%   'dataonly'     : logical, NYI (reserved for returning data without plotting).
%   'labelcones'   : logical, label cones by type in activation map (default false).
%   'lens'         : Lens object; if empty a default human lens is created
%                    and matched to cmosaic.wave for spectral plots.
%   % Line-cut specific:
%   'xdeg'         : numeric, x-position (deg) for vertical line cut (default 0).
%   'ydeg'         : numeric, y-position (deg) for horizontal line cut (default 0).
%   'thickness'    : numeric, line ROI thickness in deg (default 0.1).
%
% Outputs
%   uData  : struct with fields depending on plotType. Common fields:
%            .allE                (for activation maps/overlays)
%            .roi, .roiE, .roiIdx (for ROI/line plots)
%            .lineDistance        (for line ROI distance coordinate)
%            For 'spectralqe'     : spectral QE matrix (wave-by-coneType)
%   hdl    : matlab.ui.Figure handle of the plot. Use get(hdl,'CurrentAxes')
%            to access axes.
%
% Plot types (plotType)
%   'mosaic' | 'conemosaic'   : Render the cMosaic (delegates to cmosaic.visualize).
%
%   'excitations' | 'activations'
%       Image of activations with vertical colorbar; supports cone labels.
%       Uses fields in params passed to cmosaic.visualize().
%       Name-Values honored: 'plottitle','hdl','labelcones'
%
%   'excitations horizontal line'
%       Horizontal line ROI across mosaic at y = 'ydeg' (deg). Plots per-cone
%       excitations along the line, optionally filtered by 'conetype'.
%       Name-Values honored: 'ydeg','thickness','roi','conetype'
%
%   'excitations vertical line'
%       Vertical line ROI across mosaic at x = 'xdeg' (deg). Plots per-cone
%       excitations along the line, optionally filtered by 'conetype'.
%       Name-Values honored: 'xdeg','thickness','roi','conetype'
%
%   'excitations roi'
%       Plot excitations distribution inside a user-supplied ROI. For
%       'ellipse' ROI: histograms by cone type. For 'line' ROI: scatter of
%       excitations vs distance from the ROI origin. (Rect NYI.)
%       Name-Values honored: 'roi','conetype'
%
%   'roi'
%       Overlay the ROI outline on top of the activation image.
%       Name-Values honored: 'roi','plottitle'
%
%   'spectralqe'
%       Plot wavelength-dependent spectral quantum efficiency including lens
%       transmittance (uses provided 'lens' or creates default).
%
%   'pigmentquantalefficiency'
%       Plot photopigment quantal efficiency vs wavelength.
%
% Notes
%   • This function uses ieParamFormat(plotType) for flexible matching.
%   • Some plots expect 'allE' to be precomputed activations with time
%     matched to cmosaic.integrationTime.
%   • Dependencies: regionOfInterest, ieNewGraphWin, Lens, cmosaic.visualize,
%     cmosaic.excitations, cmosaic.coneRFpositionsDegs, cmosaic.qe.
%
% Examples
%   % Activation image
%   cmosaic.plot('excitations', allE, 'plottitle','Activation map');
%
%   % Activiation image with color labeled cones
%   cmosaic.plot('excitations', allE, 'label cones',true,'plottitle','Activation map');
%
%   % Horizontal line cut through y = 0 deg for L/M/S
%   cmosaic.plot('excitations horizontal line', allE, 'ydeg', 0);
%
%   % Vertical line cut at x = 0.1 deg for only L cones
%   cmosaic.plot('excitations vertical line', allE, 'xdeg', 0.1, 'conetype','l');
%
%   % ROI overlay and ROI histogram (ellipse ROI)
%   cmosaic.plot('roi', allE, 'roi', myEllipseROI);
%   cmosaic.plot('excitations roi', allE, 'roi', myEllipseROI);
%
%   % Spectral QE with default lens
%   cmosaic.plot('spectralqe', []);
%
% See also
%   cmosaic.visualize, cmosaic.excitations, regionOfInterest, Lens,
%   t_cMosaicBasic

%% Input parser

if ~exist('allE','var'), allE = []; end

varargin = ieParamFormat(varargin);

p = inputParser;
p.addRequired('cmosaic',@(x)(isa(x,'cMosaic')));
p.addRequired('plotType',@ischar);
p.addRequired('allE',@isnumeric);

% Excitations if precomputed
p.addParameter('conetype',{'l','m','s'},@(x)(ischar(x) || iscell(x)));
p.addParameter('roi',[],@(x)(isa(x,'regionOfInterest')));
p.addParameter('plottitle','',@ischar);
p.addParameter('dataonly',false,@islogical);   % NYI
p.addParameter('labelcones',false,@islogical)  % for excitations, show cones colored by type
p.addParameter('lens',[],@(x)(isa(x,'Lens')));
p.addParameter('coneslinewidth',2,@isnumeric);
p.addParameter('domain', 'degrees', @(x)(ischar(x) && (ismember(x, {'degrees', 'microns'}))));

% Horizontal line key val pairs
p.addParameter('xdeg',0,@isnumeric);
p.addParameter('ydeg',0,@isnumeric);
p.addParameter('thickness',0.1,@isnumeric);
p.addParameter('hdl',[],@(x)(isa(x,'matlab.ui.Figure') || isempty(x)));

p.parse(cmosaic,plotType,allE,varargin{:});

hdl    = p.Results.hdl;
pTitle = p.Results.plottitle;

% Force cone type to a cell array
conetype = p.Results.conetype;
if ischar(conetype), conetype = {conetype}; 
end

%% Different types of plots

switch ieParamFormat(plotType)
    case {'conemosaic','mosaic'}
        % We should enable passing in params
        tmp = cmosaic.visualize;        
        hdl = tmp.axesHandle;
        frame = getframe(hdl);
        uData.img = frame.cdata;
        uData.xlim = get(hdl,'xlim');
        uData.ylim = get(hdl,'ylim');

    case {'excitations','activations'}
        % Show the activations in an image
        % We should choose one - excitations or activations - for
        % consistency. 

        % Maybe we want to select out by conetype also?
        
        params = cmosaic.visualize('params');
        params.activation = allE;
        params.plotTitle = pTitle;
        params.verticalActivationColorBar = true;
        params.figureHandle = hdl;
        params.coneslinewidth = p.Results.coneslinewidth;  % Thickness of surrounding lines
        params.labelConesInActivationMap = p.Results.labelcones;
        params.domain = p.Results.domain;

        % Return
        tmp = cmosaic.visualize(params);
        hdl = tmp.figureHandle;

        % We might set a figure size (normalized position) here.
        ieFigureResize(hdl);
        uData.allE = allE;
        
    case 'excitationshorizontalline'
        % Plot the excitations along an entire horizontal line in the
        % excitations
        %
        % You can select a single cone type.  Otherwise all three
        % are plotted in RGB colors.
        xminDeg = cmosaic.eccentricityDegs(1) - cmosaic.sizeDegs(1)/2;
        xmaxDeg = cmosaic.eccentricityDegs(1) + cmosaic.sizeDegs(1)/2;
        yDeg = p.Results.ydeg;
        
        if isempty(p.Results.roi)
            roi = regionOfInterest('shape', 'line', ...
                'from', [xminDeg,yDeg], 'to', [xmaxDeg,yDeg], ...
                'thickness', p.Results.thickness);
        else
            roi = p.Results.roi;
        end
        
        hdl = ieNewGraphWin;
        roiE = cell(numel(conetype),1);
        roiIdx = cell(numel(conetype),1);
        for ii = 1:numel(conetype)
            [roiE{ii}, roiIdx{ii}] = cmosaic.excitations('roi',roi,...
                'conetype',conetype{ii},...
                'all excitations',allE);
            
            % The positions of the cones in the ROI
            pos = cmosaic.coneRFpositionsDegs(roiIdx{ii},:);
            hold on;
            thisP = plot(pos(:,1),squeeze(roiE{ii}),[coneColor(conetype{ii}),'o']);
            set(thisP,'MarkerFaceColor',coneColor(conetype{ii}));
        end
        
        hold off; grid on
        str = sprintf('Excitations (%.1f ms)',cmosaic.integrationTime*1e3);
        xlabel('Horizontal position (deg)'); ylabel(str); 
        if ~isempty(pTitle), title(pTitle); end
        
        % See how to get pos and roiE from above.
        for ii=1:numel(conetype)
            uData.pos{ii} = cmosaic.coneRFpositionsDegs(roiIdx{ii});
        end
        uData.roi = roi;
        uData.roiE = roiE;
        uData.roiIdx = roiIdx;
        
    case {'excitationsverticalline'}
        % Plot the excitations along an entire vertical line in the
        % excitations
        %
        % You can select a single cone type.  Otherwise all three
        % are plotted in RGB colors.
        yminDeg = cmosaic.eccentricityDegs(2) - cmosaic.sizeDegs(2)/2;
        ymaxDeg = cmosaic.eccentricityDegs(2) + cmosaic.sizeDegs(2)/2;
        xDeg = p.Results.xdeg;
        
        if isempty(p.Results.roi)
            roi = regionOfInterest('shape', 'line', ...
                'from', [xDeg,yminDeg], 'to', [xDeg,ymaxDeg], ...
                'thickness', p.Results.thickness);
        else
            roi = p.Results.roi;
        end
        
        hdl = ieNewGraphWin;
        roiE = cell(numel(conetype),1);
        roiIdx = cell(numel(conetype),1);
        for ii = 1:numel(conetype)
            [roiE{ii}, roiIdx{ii}] = cmosaic.excitations('roi',roi,...
                'conetype',conetype{ii},...
                'all excitations',allE);
            
            % The positions of the cones in the ROI
            pos = cmosaic.coneRFpositionsDegs(roiIdx{ii},:);
            hold on;
            thisP = plot(pos(:,2),squeeze(roiE{ii}),[coneColor(conetype{ii}),'o']);
            set(thisP,'MarkerFaceColor',coneColor(conetype{ii}));
        end
        
        hold off; grid on
        str = sprintf('Excitations (%.1f ms)',cmosaic.integrationTime*1e3);
        xlabel('Vertical position (deg)'); ylabel(str); 
        if ~isempty(pTitle), title(pTitle); end

        % See how to get pos and roiE from above.
        for ii=1:3
            uData.pos{ii} = cmosaic.coneRFpositionsDegs(roiIdx{ii});
        end
        uData.roi = roi;
        uData.roiE = roiE;
        uData.roiIdx = roiIdx;
        
    case {'roi'}
        % Draw the ROI superimposed on the excitations
        fillColor = [0.2 0.3 1];

        % Excitations
        [~,hdl] = cmosaic.plot('excitations',allE);
        axesHandle = get(hdl,'CurrentAxes');
        
        roi = p.Results.roi;
        % Show the ROI as a horizontal line on the activations
        roiOutline = roi.outline();
        
        % Draw it but keeping the activations image
        hold(axesHandle, 'on');
        
        patch(axesHandle, roiOutline.x, roiOutline.y, ...
            fillColor, 'FaceAlpha', 0.5, 'EdgeColor', fillColor*0.5, 'LineWidth', 1.0);                        
        
        uData.excitations = allE;
        uData.roi = roi;
    case {'excitationsroi'}
        % Make a plot of the excitations in an ROI
        roi = p.Results.roi;
        switch roi.shape
            case 'rect'
                disp('Rect NYI')
            case 'ellipse'
                hdl = ieNewGraphWin;
                for ii=1:numel(conetype)
                    roiE = cmosaic.excitations('roi',roi,...
                        'conetype',conetype{ii},...
                        'all excitations',allE);
                    histogram(roiE,'FaceColor',coneColor(conetype{ii}),...
                        'EdgeColor',coneColor(conetype{ii}), ...
                        'NumBins',20);
                    hold on;
                end
                xlabel('Excitations'); ylabel('Number of cones'); 
                grid on;
                
            case 'line'
                hdl = ieNewGraphWin;

                % Plot excitations along the line measured from 'from'                
                for ii = 1:numel(conetype)
                    [roiE, roiIdx] = cmosaic.excitations('roi',roi,...
                        'conetype',conetype{ii},...
                        'all excitations',allE);
                    pos = cmosaic.coneRFpositionsDegs(roiIdx,:);
                    linePosition = pos - roi.from;
                    lineDistance = sqrt(linePosition(:,1).^2 + linePosition(:,2).^2);
                    
                    % Set the color and collapse data w.r.t. the x-axis
                    hold on;
                    plot(lineDistance,squeeze(roiE),[coneColor(conetype{ii}),'o']);
                    
                    uData.lineDistance{ii} = lineDistance;
                    uData.roiE{ii} = roiE;
                end
                xlabel('Line position (deg - from)');
                ylabel(sprintf('Excitations %.1f',cmosaic.integrationTime));
                grid on;

            otherwise
                error('Unknown roi shape %s\n',roi.shape);
        end
        
    case {'spectralqe'}
        % The cMosaic does not ordinarily have a lens.  If the user
        % does not send in a lens, we use the default human lens.
        if isempty(hdl), hdl = ieNewGraphWin;
        else, figure(hdl);
        end

        if isempty(p.Results.lens)
            thisLens = Lens('wave',cmosaic.wave);
        else
            thisLens = p.Results.lens;
            thisLens.wave = cmosaic.wave;
        end
        lensT = thisLens.transmittance;
        
        % The qe incorporates the macular pigment density
        uData = diag(lensT)*cmosaic.qe;
        plot(cmosaic.wave,uData,'LineWidth',2);
        xlabel('Wavelength (nm)'); ylabel('Spectral quantum efficiency');
        if ~isempty(pTitle), title(pTitle); end

        grid on;
        
    case {'pigmentquantalefficiency'}
        hdl = ieNewGraphWin;
        plot(cmosaic.wave,cmosaic.pigment.quantalEfficiency,'LineWidth',2);
        xlabel('Wavelength (nm)'); ylabel('Quantal efficiency');
        if ~isempty(pTitle), title(pTitle); end

    otherwise
        error('Unknown plot type %s\n',plotType);
end


end

function str = coneColor(conetype)
% Make this a routine
switch lower(conetype)
    case 'l'
        str = 'r';
    case 'm'
        str = 'g';
    case 's'
        str = 'b';
end
end
