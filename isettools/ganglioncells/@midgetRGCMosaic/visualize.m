function [figureHandle, axesHandle] = visualize(obj, varargin)

    % Parse input
    p = inputParser;
    p.addParameter('figureHandle', [], @(x)(isempty(x)||isa(x, 'handle')));
    p.addParameter('axesHandle', [], @(x)(isempty(x)||isa(x, 'handle')));
    p.addParameter('maxVisualizedRFs', 18, @isscalar);
    p.addParameter('withSuperimposedOpticalImage', [], @(x)(isempty(x) || isstruct(x)));
    p.addParameter('withSuperimposedPSF', [], @(x)(isempty(x) || isstruct(x)));
    p.addParameter('xRangeDegs', 0.3, @(x)(isempty(x)||(numel(x)==1)));
    p.addParameter('yRangeDegs', 0.3, @(x)(isempty(x)||(numel(x)==1)));
    p.addParameter('xLimsDegs', [], @(x)(isempty(x)||(numel(x)==2)));
    p.addParameter('yLimsDegs', [], @(x)(isempty(x)||(numel(x)==2)));
    p.addParameter('fontSize', 16, @isscalar);
    p.addParameter('retinalMeridianAxesLabeling', true, @islogical);
    p.addParameter('plotTitle', '', @(x)(isempty(x) || ischar(x) || islogical(x)));
    p.parse(varargin{:});

    figureHandle = p.Results.figureHandle;
    axesHandle = p.Results.axesHandle;
    maxVisualizedRFs = p.Results.maxVisualizedRFs;
    superimposedOpticalImage = p.Results.withSuperimposedOpticalImage;
    superimposedPSF = p.Results.withSuperimposedPSF;
    xLimsDegs  = p.Results.xLimsDegs; 
    yLimsDegs  = p.Results.yLimsDegs;
    xRangeDegs = p.Results.xRangeDegs;
    yRangeDegs = p.Results.yRangeDegs;
    fontSize = p.Results.fontSize;
    plotTitle = p.Results.plotTitle;
    retinalMeridianAxesLabeling = p.Results.retinalMeridianAxesLabeling;

    % Set figure size
    if (isempty(figureHandle))
        figureHandle = figure(); clf;
        set(figureHandle, 'Position', [10 10 1250 1280], 'Color', [1 1 1]);
        axesHandle = subplot('Position', [0.04 0.05 0.95 0.93]);
    else
        if (isempty(axesHandle))
            figure(figureHandle);
            clf;
            set(figureHandle, 'Position', [10 10 1280 1280], 'Color', [1 1 1]);
            axesHandle = subplot('Position', [0.04 0.05 0.95 0.93]);
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
    theRetinalRFcenterMaps = obj.computeRetinalRFcenterMaps(...
        marginDegs, spatialSupportSamplesNum, ...
        'forRGCindices', sortedRGCindices(1:min([maxVisualizedRFs numel(sortedRGCindices)])));

    % Plot part of the input cone mosaic
    if (isempty(xLimsDegs))
        xLimsDegs = mRGCmosaicCenterDegs(1) + 0.5*xRangeDegs*[-1 1];
    end
    if (isempty(yLimsDegs))
        yLimsDegs = mRGCmosaicCenterDegs(2) + 0.5*yRangeDegs*[-1 1];
    end
    
    xyLimsDegs = min([xLimsDegs yLimsDegs]);
    if (xyLimsDegs < 0.5)
        xyTicksDegs = 0.1;
    elseif (xyLimsDegs < 1.0)
        xyTicksDegs = 0.2;
    elseif (xyLimsDegs < 2.5)
        xyTicksDegs = 0.5;
    elseif (xyLimsDegs < 5.0)
        xyTicksDegs = 1.0;
    elseif (xyLimsDegs < 10)
        xyTicksDegs = 2.0;
    else
        xyTicksDegs = 5.0;
    end


    xTicks = sign(mRGCmosaicCenterDegs(1)) * round(abs(mRGCmosaicCenterDegs(1)*10))/10 + xyTicksDegs*(-10:1:10);
    yTicks = sign(mRGCmosaicCenterDegs(2)) * round(abs(mRGCmosaicCenterDegs(2)*10))/10 + xyTicksDegs*(-10:1:10);

    obj.inputConeMosaic.visualize(...
            'figureHandle', figureHandle, ...
            'axesHandle', axesHandle, ...
            'visualizedConeAperture', 'lightCollectingArea4sigma', ...
            'visualizedConeApertureThetaSamples', 20, ...
            'withSuperimposedOpticalImage', superimposedOpticalImage, ...
            'withSuperimposedPSF', superimposedPSF, ...
            'domain', 'degrees', ...
            'domainVisualizationLimits', [xLimsDegs(1) xLimsDegs(2) yLimsDegs(1) yLimsDegs(2)], ...
            'domainVisualizationTicks', struct('x', xTicks, 'y', yTicks), ...
            'backgroundColor', [1 1 1], ...
            'fontSize', fontSize, ...
            'plotTitle', plotTitle);

    if (retinalMeridianAxesLabeling)
        if (obj.eccentricityDegs(1) ~= 0)
            % Change the x-label to display the horizontal retinal meridian 
            xlabel(axesHandle, sprintf('%s (degs)', strrep(obj.horizontalRetinalMeridian, 'meridian', 'retina')));
        else
            xlabel(axesHandle, '\leftarrow temporal retina               (degs)                   nasal retina \rightarrow ');
        end

        if (obj.eccentricityDegs(2) ~= 0)
            % Change the y-label to display the vertical retinal meridian 
            ylabel(axesHandle, sprintf('%s (degs)', strrep(obj.verticalRetinalMeridian, 'meridian', 'retina')));
        else
            ylabel(axesHandle, '\leftarrow inferior retina                         (degs)                          superior retina\rightarrow  ');
        end
    end

    hold(axesHandle, 'on');

    % Plot the RF maps
    for iRGC = 1:numel(sortedRGCindices)

        if (iRGC > maxVisualizedRFs)
            continue;
        end

        if (isempty(plotTitle))
           title(axesHandle, sprintf('visualized RFs: %d of %d', iRGC, numel(sortedRGCindices)));
        end

        % Retrieve the RGCindex
        targetRGCindex  = sortedRGCindices(iRGC);

        % Retrieve the computed retinal center RF map
        s = theRetinalRFcenterMaps{iRGC};
        theRF = s.centerRF;

        fprintf('Fitting ellipsoid to RF %d of %d. Please wait ...\n', iRGC, maxVisualizedRFs);

        % Fitting the discrete RF center cone map
        fitTheDiscreteRFcenterMap = true;

        if (fitTheDiscreteRFcenterMap)
            % Fit the discrete center RF map with an ellipsoidal Gaussian
            theFittedGaussian = RetinaToVisualFieldTransformer.fitScatterGaussianEllipsoid(...
                s.spatialSupportDegsX, s.spatialSupportDegsY, theRF,...
                s.inputConeWeights, obj.inputConeMosaic.coneRFpositionsDegs(s.inputConeIndices,:), ...
                'flatTopGaussian', ~true, ...
                'forcedOrientationDegs', [], ...
                'rangeForEllipseRcYRcXratio', [1/2 2], ...
                'forcedCentroidXYpos', obj.rgcRFpositionsDegs(targetRGCindex,:), ...
                'globalSearch', true, ...
                'multiStartsNum', 8);
        else
            % Fit the continuous center RF map with an ellipsoidal Gaussian
            theFittedGaussian = RetinaToVisualFieldTransformer.fitGaussianEllipsoid(...
                    s.spatialSupportDegsX, s.spatialSupportDegsY, theRF, ...
                    'flatTopGaussian', ~true, ...
                    'forcedOrientationDegs', [], ...
                    'rangeForEllipseRcYRcXratio', [1/1.4 1.4], ...
                    'forcedCentroidXYpos', obj.rgcRFpositionsDegs(targetRGCindex,:), ...
                    'globalSearch', true, ...
                    'multiStartsNum', 4);
        end


        % Plot the connection weights
        if (numel(s.inputConeIndices)>1)
            inputsCentroid = obj.rgcRFpositionsDegs(targetRGCindex,:);
            theCenterConnectedCones = s.inputConeIndices(find(s.inputConeWeights>0.001));
            for iInput = 1:numel(theCenterConnectedCones)
                switch (obj.inputConeMosaic.coneTypes(theCenterConnectedCones(iInput)))
                    case cMosaic.LCONE_ID
                        coneColor = obj.inputConeMosaic.lConeColor;
                    case cMosaic.MCONE_ID
                        coneColor = obj.inputConeMosaic.mConeColor;
                    case cMosaic.SCONE_ID
                        coneColor = obj.inputConeMosaic.sConeColor;
                end

                plot(axesHandle,[inputsCentroid(1) obj.inputConeMosaic.coneRFpositionsDegs(theCenterConnectedCones(iInput),1)], ...
                    [inputsCentroid(2) obj.inputConeMosaic.coneRFpositionsDegs(theCenterConnectedCones(iInput),2)], ...
                    'k-', 'LineWidth', 3.0, 'Color', coneColor*0.5);

                plot(axesHandle,[inputsCentroid(1) obj.inputConeMosaic.coneRFpositionsDegs(theCenterConnectedCones(iInput),1)], ...
                    [inputsCentroid(2) obj.inputConeMosaic.coneRFpositionsDegs(theCenterConnectedCones(iInput),2)], ...
                    'k-', 'LineWidth', 1.5, 'Color', coneColor);
                drawnow
            end
        end

        % Extract the fitted ellipsoidal Gaussian
        fittedEllipsoidMap = theFittedGaussian.ellipsoidMap;

        % Contour of the fitted ellipsoid at 1 Rc
        oneRcLevel = [1 exp(-1)];

        if (numel(theCenterConnectedCones) == 1)
            % Single input. Go down to 10%
            oneRcLevel = [1 0.1];
        end

        cMap = brewermap(256, '*oranges');
        alpha = 0.15;
        contourLineColor = [0.2 0.2 0.2];
        
        cMosaic.semiTransparentContourPlot(axesHandle, s.spatialSupportDegsX, s.spatialSupportDegsY, fittedEllipsoidMap/max(fittedEllipsoidMap(:)), ...
            oneRcLevel, cMap, alpha, contourLineColor, ...
            'lineWidth', 3.0, ...
            'edgeAlpha', 0.7);
        drawnow;

    end


end

