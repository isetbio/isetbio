function visualizeSingleRetinalRF(obj,theRGCindex, varargin)
    % Parse input
    p = inputParser;
    p.addParameter('figureHandle', [], @(x)(isempty(x)||isa(x, 'handle')));
    p.addParameter('axesHandle', [], @(x)(isempty(x)||isa(x, 'handle')));
    p.addParameter('fontSize', 16, @isscalar);
    p.addParameter('plotTitle', '', @ischar);
    p.addParameter('plotTitleColor', [0 0 0], @isnumeric);
    p.addParameter('xRangeDegs', 0.1, @isscalar);
    p.addParameter('yRangeDegs', 0.1, @isscalar);
    p.addParameter('showInputConeMosaic', true, @islogical);
    p.addParameter('showConeWeights', true, @islogical);
    p.parse(varargin{:});

    hFig = p.Results.figureHandle;
    ax = p.Results.axesHandle;
    fontSize = p.Results.fontSize;
    plotTitle = p.Results.plotTitle;
    plotTitleColor = p.Results.plotTitleColor;
    xRange = p.Results.xRangeDegs;
    yRange = p.Results.yRangeDegs;

    % Find this RGC's center input cone indices and their weights
    connectivityVector = full(squeeze(obj.rgcRFcenterConeConnectivityMatrix(:, theRGCindex)));
    inputConeIndices = find(connectivityVector > 0.0001);

    xyCenter = mean(obj.inputConeMosaic.coneRFpositionsDegs(inputConeIndices,:),1);
    xyLims  =[ xyCenter(1)-xRange*0.5 xyCenter(1)+xRange*0.5  xyCenter(2)-yRange*0.5  xyCenter(2)+yRange*0.5];

    obj.inputConeMosaic.visualize(...
            'figureHandle', hFig, ...
            'axesHandle', ax, ...
            'domain', 'degrees', ...
            'outlinedConesWithIndices', inputConeIndices, ...
            'visualizedConeAperture', 'lightCollectingAreaCharacteristicDiameter', ...
            'visualizedConeApertureThetaSamples', 30, ...
            'backgroundColor', [1 1 1], ...
            'domainVisualizationLimits', xyLims, ...
            'plotTitle', plotTitle, ...
            'plotTitleColor', plotTitleColor, ...
            'noXLabel', true, ...
            'noYLabel', true, ...
            'fontSize', fontSize);
  
    hold(ax, 'on');
    drawnow;

end

function plotConesInTheNeighborhood(ax, inputConeMosaic, xLims, yLims)

    coneXpositionDegs = inputConeMosaic.coneRFpositionsDegs(:,1);
    coneYpositionDegs = inputConeMosaic.coneRFpositionsDegs(:,2);
    coneIndicesInNeigborhood = find(...
        (coneXpositionDegs>=xLims(1)) & ...
        (coneXpositionDegs<=xLims(2)) & ...
        (coneYpositionDegs>=yLims(1)) & ...
        (coneYpositionDegs<=yLims(2)) ...
        );

    % Plot all cones in the neighborhood
    for iCone = 1:numel(coneIndicesInNeigborhood)
        coneIndex = coneIndicesInNeigborhood(iCone);
        xo = inputConeMosaic.coneRFpositionsDegs(coneIndex,1);
        yo = inputConeMosaic.coneRFpositionsDegs(coneIndex,2);
        r = inputConeMosaic.coneApertureDiametersDegs(coneIndex) * ...
            inputConeMosaic.coneApertureToConeCharacteristicRadiusConversionFactor;
        xx = xo + r * cosd(0:10:360);
        yy = yo + r * sind(0:10:360);
        switch (inputConeMosaic.coneTypes(coneIndex))
            case cMosaic.LCONE_ID
                faceColor = [1 0.2 0.5];
            case cMosaic.MCONE_ID
                faceColor = [0.2 1 0.5];
            case cMosaic.SCONE_ID
                faceColor = [0.5 0.1 0.9];
        end

        patch(xx, yy, [0 0 0], 'EdgeColor', 'none', 'FaceAlpha', 0.8, ...
            'FaceColor', faceColor, 'LineWidth', 0.5, ...
            'LineStyle', '-', ...
            'Parent', ax);
    end
end
