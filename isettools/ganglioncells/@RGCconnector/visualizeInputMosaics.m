function [hFig, ax, XLims, YLims] = visualizeInputMosaics(obj, varargin)
% Outputs:
%    hFig               - The figure handle
%    ax                 - The axes handle
%    XLims, YLims       - The limits [min max] of the x- and y-axes
%
% Optional key/value pairs
%   'figureHandle'      - The figure handle on which to render the figure
%   'axesHandle'        - The axes handle on which to render the figure
%   'inputConeMosaic'   - The input cone mosaic
%   'XLims', 'YLims'    - The limits [min max] of the x- and y-axes
%   'titleString'       - A title string
%   'thetaSamples'      - How many samples to use for rendering the RF disk
%
% History:
%   5/11/2022       NPC     Wrote it
%

    p = inputParser;
    p.addParameter('figureHandle', [], @(x)(isempty(x)||isa(x, 'handle')));
    p.addParameter('axesHandle', [], @(x)(isempty(x)||isa(x, 'handle')));
    p.addParameter('titleString', '', @(x)(isempty(x) || (ischar(x))));
    p.addParameter('thetaSamples', 20, @isnumeric);
    p.addParameter('XLims', [], @isnumeric);
    p.addParameter('YLims', [], @isnumeric);
    p.parse(varargin{:});
    
    hFig = p.Results.figureHandle;
    ax = p.Results.axesHandle;
    titleString = p.Results.titleString;
    thetaSamples = p.Results.thetaSamples;
    XLims = p.Results.XLims;
    YLims = p.Results.YLims;

    % Generate disk outline
    thetas = linspace(0,360,thetaSamples);
    diskOutline = 0.5*[cosd(thetas); sind(thetas)]';

    if (isempty(ax))
        if (isempty(hFig))
            hFig = figure(); clf;
            set(hFig, 'Color', [1 1 1], 'Position', [10 10 1700 500]);
        end
        ax = subplot('Position', [0.05 0.07 0.93 0.9]);
    end
    
    hold(ax, 'on')

    % Plot the cones
    obj.visualizeConePositions(ax, diskOutline);
  
    % Plot the RGCs
    [f,v] = RGCconnector.facesAndVertices(obj.RGCRFpositionsMicrons, obj.RGCRFspacingsMicrons, diskOutline);
    patch(ax,'Faces', f, 'Vertices', v, 'FaceColor', [0.75 0.75 0.75], 'EdgeColor', [0 0 0], ...
        'FaceAlpha', 0.35, 'LineWidth', 2.0);


    % Finalize
    axis(ax, 'equal');

    minRGCXY = min(obj.RGCRFpositionsMicrons,[],1);
    minConeXY = min(obj.inputConeMosaic.coneRFpositionsMicrons,[],1);
    minX = min([minConeXY(1) minRGCXY(1)]);
    minY = min([minConeXY(2) minRGCXY(2)]);

    maxRGCXY = max(obj.RGCRFpositionsMicrons,[],1);
    maxConeXY = max(obj.inputConeMosaic.coneRFpositionsMicrons,[],1);
    maxX = max([maxConeXY(1) maxRGCXY(1)]);
    maxY = max([maxConeXY(2) maxRGCXY(2)]);

    minX = minRGCXY(1); maxX = maxRGCXY(1);
    minY = minRGCXY(2); maxY = maxRGCXY(2);

    
    maxSpacing = 0.5*max(obj.RGCRFspacingsMicrons);
    if (isempty(XLims))
        XLims = [minX-maxSpacing maxX+maxSpacing];
    end
    if (isempty(YLims))
        YLims = [minY-maxSpacing maxY+maxSpacing];
    end

    set(ax, 'XLim', XLims, 'YLim', YLims, 'FontSize', 16);

    if (~isempty(titleString))
        title(ax, titleString);
    end

    box(ax, 'on')
    xlabel(ax, 'microns'); 
    ylabel(ax, 'microns');
    drawnow;
end