function visualize(obj, varargin)

    % Parse input
    p = inputParser;
    p.addParameter('figureHandle', [], @(x)(isempty(x)||isa(x, 'handle')));
    p.addParameter('axesHandle', [], @(x)(isempty(x)||isa(x, 'handle')));
    p.addParameter('maxVisualizedRFs', 7, @isscalar);
    p.addParameter('xLimsDegs', [], @(x)(isempty(x)||(numel(x)==2)));
    p.addParameter('yLimsDegs', [], @(x)(isempty(x)||(numel(x)==2)));
    p.addParameter('xRangeDegs', [], @(x)(isempty(x)||(numel(x)==2)));
    p.addParameter('yRangeDegs', [], @(x)(isempty(x)||(numel(x)==2)));
    p.addParameter('fontSize', 16, @isscalar);
    p.addParameter('retinalMeridianAxesLabeling', false, @islogical);
    p.addParameter('plotTitle', '', @(x)(isempty(x) || ischar(x) || islogical(x)));
    p.parse(varargin{:});

    figureHandle = p.Results.figureHandle;
    axesHandle = p.Results.axesHandle;
    maxVisualizedRFs = p.Results.maxVisualizedRFs;
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
    if (isempty(xLimsDegs))
        xLimsDegs = mRGCmosaicCenterDegs(1) + 0.15*[-1 1];
    end
    if (isempty(yLimsDegs))
        yLimsDegs = mRGCmosaicCenterDegs(2) + 0.15*[-1 1];
    end
    
    if (~isempty(xRangeDegs))
        xLimsDegs = mRGCmosaicCenterDegs(1) + xRangeDegs;
    end

    if (~isempty(yRangeDegs))
        yLimsDegs = mRGCmosaicCenterDegs(2) + yRangeDegs;
    end

    xTicks = sign(mRGCmosaicCenterDegs(1)) * round(abs(mRGCmosaicCenterDegs(1)*10))/10 + 0.1*(-2:1:2);
    yTicks = sign(mRGCmosaicCenterDegs(2)) * round(abs(mRGCmosaicCenterDegs(2)*10))/10 + 0.1*(-2:1:2);

    obj.inputConeMosaic.visualize(...
            'figureHandle', figureHandle, ...
            'axesHandle', axesHandle, ...
            'visualizedConeAperture', 'lightCollectingArea4sigma', ...
            'visualizedConeApertureThetaSamples', 20, ...
            'domain', 'degrees', ...
            'domainVisualizationLimits', [xLimsDegs(1) xLimsDegs(2) yLimsDegs(1) yLimsDegs(2)], ...
            'domainVisualizationTicks', struct('x', xTicks, 'y', yTicks), ...
            'backgroundColor', [1 1 1], ...
            'fontSize', fontSize, ...
            'plotTitle', plotTitle);

    if (retinalMeridianAxesLabeling)
        if (obj.eccentricityDegs(1) ~= 0)
            % Change the x-label to display the horizontal retinal meridian 
            xlabel(axesHandle, sprintf('%s', strrep(obj.horizontalRetinalMeridian, 'meridian', 'retina')));
        else
            xlabel(axesHandle, '\leftarrow nasal retina    \rightarrow temporal retina');
        end

        if (obj.eccentricityDegs(2) ~= 0)
            % Change the y-label to display the vertical retinal meridian 
            ylabel(axesHandle, sprintf('%s', strrep(obj.verticalRetinalMeridian, 'meridian', 'retina')));
        else
            ylabel(axesHandle, '\leftarrow inferior retina    \rightarrow superior retina');
        end
    end

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

