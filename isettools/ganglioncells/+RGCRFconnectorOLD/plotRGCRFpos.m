function [hFig, ax, XLims, YLims] = plotRGCRFpos(RGCRFposMicrons, varargin)
% Plot the RF positions of the RGC mosaic, and optionally of the input cone mosaic
%
% Syntax:
%   RGCRFconnector.plotRGCRFpos(RGCRFposMicrons, varargin);
%
% Description:
%   Plot the RF positions of the RGC mosaic and, optionally of the input 
%   cone mosaic
%
% Inputs:
%    RGCRFposMicrons    - The positions of the RGC RFs in microns
%
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
    p.addParameter('inputConeMosaic', [], @(x)(isempty(x)||(isa(x, 'cMosaic'))));
    p.addParameter('titleString', '', @(x)(isempty(x) || (ischar(x))));
    p.addParameter('thetaSamples', 15, @isnumeric);
    p.addParameter('XLims', [], @isnumeric);
    p.addParameter('YLims', [], @isnumeric);
    p.parse(varargin{:});
    
    hFig = p.Results.figureHandle;
    ax = p.Results.axesHandle;
    inputConeMosaic = p.Results.inputConeMosaic;
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
            set(hFig, 'Color', [1 1 1]);
        end
        ax = subplot('Position', [0.1 0.1 0.8 0.8]);
    end
    
    hold(ax, 'on')
    if (~isempty(inputConeMosaic))
        % Plot the cones
        coneTypes = [cMosaic.LCONE_ID cMosaic.MCONE_ID cMosaic.SCONE_ID];
        coneColors = [inputConeMosaic.lConeColor; inputConeMosaic.mConeColor; inputConeMosaic.sConeColor];
        
        for iConeType = 1:numel(coneTypes)
            idx = find(inputConeMosaic.coneTypes == coneTypes(iConeType));
            allConeRFpositionsMicrons = inputConeMosaic.coneRFpositionsMicrons(idx,:);
            allConeRFspacingsMicrons = inputConeMosaic.coneRFspacingsMicrons(idx);
            [f,v] = facesAndVertices(allConeRFpositionsMicrons, allConeRFspacingsMicrons, diskOutline);
            theColor = squeeze(coneColors(iConeType,:));
            patch(ax,'Faces', f, 'Vertices', v, 'FaceColor', theColor, 'EdgeColor', theColor*0.5);
        end
    end

    % Plot the RGCRFs
    RGCRFspacingsMicrons = RGCmodels.Watson.convert.positionsToSpacings(RGCRFposMicrons);
    [f,v] = facesAndVertices(RGCRFposMicrons, RGCRFspacingsMicrons, diskOutline);
    patch(ax,'Faces', f, 'Vertices', v, 'FaceColor', [0.5 0.5 0.5], 'EdgeColor', [0.3 0.3 0.3], 'FaceAlpha', 0.35, 'LineWidth', 1.0);

    % Finalize
    axis(ax, 'equal');

    if (isempty(XLims))
        XLims = [min(RGCRFposMicrons(:,1))-max(RGCRFspacingsMicrons) max(RGCRFposMicrons(:,1))+max(RGCRFspacingsMicrons)];
    end
    if (isempty(YLims))
        YLims = [min(RGCRFposMicrons(:,2))-max(RGCRFspacingsMicrons) max(RGCRFposMicrons(:,2))+max(RGCRFspacingsMicrons)];
    end

    set(ax, 'XLim', XLims, 'YLim', YLims, 'FontSize', 16);

    if (~isempty(titleString))
        title(ax, titleString);
    end

end

function [f,v] = facesAndVertices(positions, spacings, diskOutline)
    thetaSamples = size(diskOutline,1);
    rfsNum = size(positions, 1);
    X = zeros(rfsNum*thetaSamples,1);
    Y = X;
    f = zeros(rfsNum,thetaSamples);
    for iRF = 1:rfsNum
        ii = (iRF-1)*thetaSamples + (1:thetaSamples);
        f(iRF,:) = ii;
        X(ii,:) = positions(iRF,1) + spacings(iRF) * diskOutline(:,1);
        Y(ii,:) = positions(iRF,2) + spacings(iRF) * diskOutline(:,2);
    end
    v = [X(:) Y(:)];
end
