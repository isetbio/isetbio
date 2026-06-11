function [hFig, ax, XLims, YLims] = visualizeInputLattices(obj, varargin)
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

    thetas = linspace(0,360,4);
    triangleOutline = 0.5*[cosd(thetas); sind(thetas)]';

    thetas60 = linspace(0,360,17);
    diskOutline = 0.5*[cosd(thetas60); sind(thetas60)]';

    if (isempty(ax))
        if (isempty(hFig))
            hFig = figure(); clf;
            set(hFig, 'Color', [1 1 1], 'Position', [10 10 800 500]);
        end
        ax = subplot('Position', [0.05 0.07 0.93 0.88]);
    end
    
    hold(ax, 'on')

    % Plot source lattice positions in black triangles
    [f,v] = MosaicConnector.facesAndVertices(...
        obj.sourceLattice.RFpositionsMicrons, ...
        obj.sourceLattice.RFspacingsMicrons, triangleOutline);
    patch(ax,'Faces', f, 'Vertices', v, 'FaceColor', [0.25 0.25 0.25], 'EdgeColor', [0 0 0], ...
        'FaceAlpha', 0.7, 'LineWidth', 2.0, 'LineStyle', '-');


    if (isempty(obj.destinationRFspacingsFromCentroids))
        % Plot original destination lattice positions in cyan disks
        [f,v] = MosaicConnector.facesAndVertices(...
            obj.destinationLattice.RFpositionsMicrons, ...
            obj.destinationLattice.RFspacingsMicrons, diskOutline);
        faceColor = [0.99 0.55 0.65];
        edgeColor = faceColor*0.5;
    else
        % Plot the inputRF-based destination RF positions in orange disks
        [f,v] = MosaicConnector.facesAndVertices(...
            obj.destinationRFcentroidsFromInputs, ...
            obj.destinationRFspacingsFromCentroids, diskOutline);
        faceColor = [0.25 0.85 0.99];
        edgeColor = faceColor*0.5;
    end

    patch(ax,'Faces', f, 'Vertices', v, 'FaceColor', faceColor, 'EdgeColor', edgeColor, ...
        'FaceAlpha', 0.75, 'LineWidth', 2.0, 'LineStyle', '-');


    % Finalize
    axis(ax, 'equal');

    minXYdest = min(obj.destinationLattice.RFpositionsMicrons,[],1);
    %minXYsrc =  min(obj.sourceLattice.RFpositionsMicrons,[],1);

    maxXYdest = max(obj.destinationLattice.RFpositionsMicrons,[],1);
    %maxXYsrc =  max(obj.sourceLattice.RFpositionsMicrons,[],1);

    maxDestSpacing = max(obj.destinationLattice.RFspacingsMicrons);
    if (isempty(XLims))
        XLims = [minXYdest(1) maxXYdest(1)] + maxDestSpacing*[-1 1];
    end
    if (isempty(YLims))
        YLims = [minXYdest(2) maxXYdest(2)] + maxDestSpacing*[-1 1];
    end
   
    set(ax, 'XLim', XLims, 'YLim', YLims, 'FontSize', 16);

    if (~isempty(titleString))
        title(ax, titleString);
    end

    title(sprintf('cyan hexagons: %s, gray triangles: %s', ...
        obj.destinationLattice.name, obj.sourceLattice.name));
    box(ax, 'on')
    xlabel(ax, 'microns'); 
    ylabel(ax, 'microns');
    drawnow;
end