function visualize(obj, varargin)

    % Parse input
    p = inputParser;
    p.addParameter('figureHandle', [], @(x)(isempty(x)||isa(x, 'handle')));
    p.addParameter('axesHandle', [], @(x)(isempty(x)||isa(x, 'handle')));
    p.addParameter('maxVisualizedRFs', 7, @isscalar);
    p.addParameter('xLims', [], @(x)(isempty(x)||(numel(x)==2)));
    p.addParameter('yLims', [], @(x)(isempty(x)||(numel(x)==2)));
    p.addParameter('fontSize', 16, @isscalar);
    p.parse(varargin{:});

    figureHandle = p.Results.figureHandle;
    axesHandle = p.Results.axesHandle;
    maxVisualizedRFs = p.Results.maxVisualizedRFs;
    xLims  = p.Results.xLims; 
    yLims  = p.Results.yLims;
    fontSize = p.Results.fontSize;

    % Set figure size
    if (isempty(figureHandle))
        figureHandle = figure(); clf;
        set(figureHandle, 'Position', [10 10 700 700], 'Color', [1 1 1]);
        axesHandle = subplot('Position', [0.11 0.07 0.85 0.90]);
    else
        if (isempty(axesHandle))
            figure(figureHandle);
            clf;
            set(figureHandle, 'Position', [10 10 700 700], 'Color', [1 1 1]);
            axesHandle = subplot('Position', [0.11 0.07 0.85 0.90]);
        end
        cla(axesHandle);
    end


    % Compute center of mosaic
    mRGCmosaicCenterDegs = mean(obj.rgcRFpositionsDegs,1);

    % Sort RGCs according to their distance from the mosaic center
    ecc = sum((bsxfun(@minus, obj.rgcRFpositionsDegs, mRGCmosaicCenterDegs)).^2,2);
    [~,sortedRGCindices] = sort(ecc, 'ascend');

    % Compute the retinal RFcenter maps
    marginDegs = min([0.5 0.4*min(obj.sizeDegs)]);
    spatialSupportSamplesNum = 256;
    theRetinalRFcenterMaps = obj.computeRetinalRFcenterMaps(marginDegs, spatialSupportSamplesNum);

    % Plot part of the input cone mosaic
    if (isempty(xLims))
        xLims = mRGCmosaicCenterDegs(1) + 0.15*[-1 1];
    end
    if (isempty(yLims))
        yLims = mRGCmosaicCenterDegs(2) + 0.5*[-1 1];
    end
    
    xTicks = sign(mRGCmosaicCenterDegs(1)) * round(abs(mRGCmosaicCenterDegs(1)*10))/10 + 0.1*(-2:1:2);
    yTicks = sign(mRGCmosaicCenterDegs(2)) * round(abs(mRGCmosaicCenterDegs(2)*10))/10 + 0.1*(-2:1:2);

    obj.inputConeMosaic.visualize(...
            'figureHandle', figureHandle, ...
            'axesHandle', axesHandle, ...
            'visualizedConeAperture', 'lightCollectingArea4sigma', ...
            'visualizedConeApertureThetaSamples', 20, ...
            'domain', 'degrees', ...
            'domainVisualizationLimits', [xLims(1) xLims(2) yLims(1) yLims(2)], ...
            'domainVisualizationTicks', struct('x', xTicks, 'y', yTicks), ...
            'backgroundColor', [1 1 1], ...
            'fontSize', fontSize);

    hold(axesHandle, 'on');

    % Plot the RF maps
    for iRGC = 1:numel(sortedRGCindices)

        if (iRGC > maxVisualizedRFs)
            continue;
        end

        % Retrieve the RGCindex
        targetRGCindex  = sortedRGCindices(iRGC);

        % Retrieve the computed retinal center RF map
        s = theRetinalRFcenterMaps{targetRGCindex};
        theRF = s.centerRF;

        % Fit the continuous center RF map with an ellipsoidal Gaussian
        theFittedGaussian = RetinaToVisualFieldTransformer.fitGaussianEllipsoid(...
                s.spatialSupportDegsX, s.spatialSupportDegsY, theRF, ...
                'flatTopGaussian', ~true, ...
                'forcedOrientationDegs', [], ...
                'rangeForEllipseRcYRcXratio', [1/1.4 1.4], ...
                'forcedCentroidXYpos', obj.rgcRFpositionsDegs(targetRGCindex,:), ...
                'globalSearch', true, ...
                'multiStartsNum', 8);


        % Plot the connection weights
        if (numel(s.inputConeIndices)>1)
            inputsCentroid = obj.rgcRFpositionsDegs(targetRGCindex,:);
            for iInput = 1:numel(s.inputConeIndices)
                if (s.inputConeWeights(iInput)>0.01)
                    switch (obj.inputConeMosaic.coneTypes(s.inputConeIndices(iInput)))
                        case cMosaic.LCONE_ID
                            coneColor = obj.inputConeMosaic.lConeColor;
                        case cMosaic.MCONE_ID
                            coneColor = obj.inputConeMosaic.mConeColor;
                        case cMosaic.SCONE_ID
                            coneColor = obj.inputConeMosaic.sConeColor;
                    end

                    plot(axesHandle,[inputsCentroid(1) obj.inputConeMosaic.coneRFpositionsDegs(s.inputConeIndices(iInput),1)], ...
                            [inputsCentroid(2) obj.inputConeMosaic.coneRFpositionsDegs(s.inputConeIndices(iInput),2)], ...
                            'k-', 'LineWidth', 3.0, 'Color', coneColor*0.5);

                    plot(axesHandle,[inputsCentroid(1) obj.inputConeMosaic.coneRFpositionsDegs(s.inputConeIndices(iInput),1)], ...
                            [inputsCentroid(2) obj.inputConeMosaic.coneRFpositionsDegs(s.inputConeIndices(iInput),2)], ...
                            'k-', 'LineWidth', 1.5, 'Color', coneColor);
                    drawnow
                end
            end
        end

        % Extract the fitted ellipsoidal Gaussian
        fittedEllipsoidMap = theFittedGaussian.ellipsoidMap;

        % Contour of the fitted ellipsoid at 1 Rc
        oneRcLevel = [1 exp(-1)];
        cMap = brewermap(256, '*oranges');
        alpha = 0.25;
        contourLineColor = [0.2 0.2 0.2];
        
        cMosaic.semiTransparentContourPlot(axesHandle, s.spatialSupportDegsX, s.spatialSupportDegsY, fittedEllipsoidMap/max(fittedEllipsoidMap(:)), ...
            oneRcLevel, cMap, alpha, contourLineColor, ...
            'lineWidth', 3.0, ...
            'edgeAlpha', 0.7);
        drawnow;

    end


end

