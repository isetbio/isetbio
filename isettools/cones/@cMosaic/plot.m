function [uData, hdl] = plot(cmosaic,plotType, allE, varargin)
% plot methods for the cMosaic
%
% Syopsis
%    [uData, hdl] = plot(cmosaic, plotType, allE, varargin)
%
% Inputs
%    cmosaic   - cMosaic class
%    plotType  - See below.  Many.
%    
% Optional key/val pairs
%    roi 
%    cone type
%    hdl - Figure handle (matlab.ui.Figure).  Used for overlaying
%           curves.
%
% Output
%    uData - Struct with the plot data including the ROI
%    hdl   - Plot figure handle.  Use get(hdl,'CurrentAxes') for axis
%
% Plot types
%
%    excitations - pull out excitations of various types
%       cmosaic.plot('excitations',allE) 
%       cmosaic.plot('excitations horizontal line',allE,'ydeg',0)
%       cmosaic.plot('excitations roi',allE,'cone type',conetype)
%
%    roi - show the ROI superimposed on the excitation image
%        cmosaic.plot('roi',allE,'roi',regionOfInterest)
%
%    plot title - Logical.  Show the title or not.
%
% See also
%   t_cMosaicBasic
%

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

p.addParameter('lens',[],@(x)(isa(x,'Lens')));

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
    case {'excitations','activations'}
        % Show the activations in an image
        % We should choose one - excitations or activations - for
        % consistency. 

        % Maybe we want to select out by conetype also?
        
        params = cmosaic.visualize('params');
        params.activation = allE;
        params.plotTitle = pTitle;
        params.verticalActivationColorBar = true;
        
        % Return
        tmp = cmosaic.visualize(params);
        hdl = tmp.figureHandle;
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
        for ii=1:3
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
