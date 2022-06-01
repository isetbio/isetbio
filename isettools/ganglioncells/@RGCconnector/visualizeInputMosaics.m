function [hFig, ax, XLims, YLims] = visualizeInputMosaics(obj, varargin)
% Outputs:
%    hFig               - The figure handle
%    ax                 - The axes handle
%    XLims, YLims       - The limits [min max] of the x- and y-axes
%
% Optional key/value pairs
%   'figureHandle'      - The figure handle on which to render the figure
%   'axesHandle'        - The axes handle on which to render the figure
%   'XLims', 'YLims'    - The limits [min max] of the x- and y-axes
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
    thetaSamples = p.Results.thetaSamples;
    titleString = p.Results.titleString;
    XLims = p.Results.XLims;
    YLims = p.Results.YLims;

    % Generate disk outline
    thetas = linspace(0,360,thetaSamples);
    diskOutline = 0.5*[cosd(thetas); sind(thetas)]';
    thetas60 = linspace(0,360,7);
    hexOutline = 0.5*[cosd(thetas60); sind(thetas60)]';

    if (isempty(ax))
        if (isempty(hFig))
            hFig = figure(); clf;
            set(hFig, 'Color', [1 1 1], 'Position', [10 10 1700 500]);
        end
        ax = subplot('Position', [0.05 0.07 0.93 0.9]);
    end
    
    hold(ax, 'on')

    % Plot original RGCRF positions in gray hegagons
    [f,v] = RGCconnector.facesAndVertices(...
        obj.RGCRFpositionsMicrons, ...
        obj.RGCRFspacingsMicrons, hexOutline);
    patch(ax,'Faces', f, 'Vertices', v, 'FaceColor', [0.25 0.25 0.25], 'EdgeColor', [0 0 0], ...
        'FaceAlpha', 0.7, 'LineWidth', 2.0, 'LineStyle', '-');

    
    % Plot the centroid positions of RGCs with cone inputs in cyan disks
    idx = find(obj.RGCRFcentroidsFromInputs(:,1) ~= Inf);
    [f,v] = RGCconnector.facesAndVertices(...
        obj.RGCRFcentroidsFromInputs(idx,:), ...
        obj.RGCRFspacingsMicrons(idx), diskOutline);
    patch(ax,'Faces', f, 'Vertices', v, 'FaceColor', [0.25 0.75 0.75], 'EdgeColor', [0 0 0], ...
        'FaceAlpha', 0.7, 'LineWidth', 2.0, 'LineStyle', '-');

    idx = find(obj.RGCRFcentroidsFromInputs(:,1) == Inf);
    [f,v] = RGCconnector.facesAndVertices(...
        obj.RGCRFcentroidsFromInputs(idx,:), ...
        obj.RGCRFspacingsMicrons(idx), diskOutline);
    patch(ax,'Faces', f, 'Vertices', v, 'FaceColor', [0.25 0.25 0.25], 'EdgeColor', [0 0 0], ...
        'FaceAlpha', 0.7, 'LineWidth', 2.0, 'LineStyle', '-');



    % Finalize
    axis(ax, 'equal');

    minRGCXY = min(obj.RGCRFpositionsMicrons,[],1);
    minConeXY = min(obj.inputConeMosaic.coneRFpositionsMicrons,[],1);

    maxRGCXY = max(obj.RGCRFpositionsMicrons,[],1);
    maxConeXY = max(obj.inputConeMosaic.coneRFpositionsMicrons,[],1);

    minX = 0.1*minConeXY(1) + 0.9*minRGCXY(1); 
    maxX = 0.1*maxConeXY(1) + 0.9* maxRGCXY(1);
    minY = 0.1*minConeXY(2) + 0.9* minRGCXY(2);
    maxY = 0.1*maxConeXY(2) + 0.9* maxRGCXY(2);
    
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

    title('hexagons: original RF lattice, disks: cone-input based RF lattice')
    box(ax, 'on')
    xlabel(ax, 'microns'); 
    ylabel(ax, 'microns');
    drawnow;
end