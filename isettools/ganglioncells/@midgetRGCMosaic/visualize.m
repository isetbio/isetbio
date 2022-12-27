function [figureHandle, axesHandle] = visualize(obj, varargin)

    % Parse input
    p = inputParser;
    p.addParameter('figureHandle', [], @(x)(isempty(x)||isa(x, 'handle')));
    p.addParameter('axesHandle', [], @(x)(isempty(x)||isa(x, 'handle')));
    p.addParameter('eccentricitySamplingGrid', [], @(x)(isempty(x) || (isnumeric(x) && (size(x,2) == 2)) ));
    p.addParameter('maxVisualizedRFs', 0, @(x)(isempty(x) || isscalar(x)));
    p.addParameter('showConnectionsToCones', true, @islogical);
    p.addParameter('withSuperimposedConeMosaic', true, @islogical);
    p.addParameter('withSuperimposedOpticalImage', [], @(x)(isempty(x) || isstruct(x)));
    p.addParameter('withSuperimposedPSF', [], @(x)(isempty(x) || isstruct(x)));
    p.addParameter('inputPoolingVisualization', '', @(x)((isempty(x))||(ismember(x, {'centerOnly', 'surroundOnly', 'centerAndSurround'}))));
    p.addParameter('xRangeDegs', [], @(x)(isempty(x)||(numel(x)==1)));
    p.addParameter('yRangeDegs', [], @(x)(isempty(x)||(numel(x)==1)));
    p.addParameter('xLimsDegs', [], @(x)(isempty(x)||isinf(x)||(numel(x)==2)));
    p.addParameter('yLimsDegs', [], @(x)(isempty(x)||isinf(x)||(numel(x)==2)));
    p.addParameter('domainVisualizationTicks', [], @(x)(isempty(x)||(isstruct(x)&&((isfield(x, 'x'))&&(isfield(x,'y'))))));
    p.addParameter('noXLabel', false, @islogical);
    p.addParameter('noYLabel', false, @islogical);
    p.addParameter('contourLineWidth', 3, @isscalar);
    p.addParameter('fontSize', 16, @isscalar);
    p.addParameter('retinalMeridianAxesLabeling', true, @islogical);
    p.addParameter('plotTitle', '', @(x)(isempty(x) || ischar(x) || islogical(x)));
    p.parse(varargin{:});

    figureHandle = p.Results.figureHandle;
    axesHandle = p.Results.axesHandle;
    eccentricitySamplingGrid = p.Results.eccentricitySamplingGrid;
    maxVisualizedRFs = p.Results.maxVisualizedRFs;
    showConnectionsToCones = p.Results.showConnectionsToCones;
    superimposedOpticalImage = p.Results.withSuperimposedOpticalImage;
    superimposedPSF = p.Results.withSuperimposedPSF;
    superimposedConeMosaic = p.Results.withSuperimposedConeMosaic;
    xRangeDegs = p.Results.xRangeDegs;
    yRangeDegs = p.Results.yRangeDegs;
    xLimsDegs  = p.Results.xLimsDegs; 
    yLimsDegs  = p.Results.yLimsDegs;
    inputPoolingVisualization = p.Results.inputPoolingVisualization;
    domainVisualizationTicks = p.Results.domainVisualizationTicks;
    contourLineWidth = p.Results.contourLineWidth;
    noXLabel = p.Results.noXLabel;
    noYLabel = p.Results.noYLabel;
    fontSize = p.Results.fontSize;
    plotTitle = p.Results.plotTitle;
    retinalMeridianAxesLabeling = p.Results.retinalMeridianAxesLabeling;

    % Set figure size
    if (isempty(figureHandle))
        figureHandle = figure(); clf;
        set(figureHandle, 'Position', [10 10 1250 1280], 'Color', [1 1 1]);
        axesHandle = subplot('Position', [0.06 0.07 0.93 0.9]);
    else
        if (isempty(axesHandle))
            figure(figureHandle);
            clf;
            set(figureHandle, 'Position', [10 10 1280 1280], 'Color', [1 1 1]);
            axesHandle = subplot('Position', [0.06 0.07 0.93 0.9]);
        end
        cla(axesHandle);
    end

    % Compute center of mosaic
    mRGCmosaicCenterDegs = mean(obj.rgcRFpositionsDegs,1);

    % Sort RGCs according to their distance from the mosaic center
    ecc = sum((bsxfun(@minus, obj.rgcRFpositionsDegs, mRGCmosaicCenterDegs)).^2,2);
    [~,sortedRGCindices] = sort(ecc, 'ascend');

    if (maxVisualizedRFs < 0)
        % If negative half of the visualizedRFs will be from the mosaic
        % center and the other half from its 4 corners
        maxVisualizedRFs = -maxVisualizedRFs;
        m1 = min([floor(maxVisualizedRFs/2) numel(sortedRGCindices)]);
        sortedRGCindices = unique([sortedRGCindices(1:m1); sortedRGCindices(end:-1:numel(sortedRGCindices)-m1)]);
    end

    % Compute the retinal RFcenter maps
    if (isempty(inputPoolingVisualization))
        marginDegs = min([0.5 0.4*min(obj.sizeDegs)]);
        if (maxVisualizedRFs > 0)
            spatialSupportSamplesNum = 256;
            theRetinalRFcenterMaps = obj.computeRetinalRFcenterMaps(...
                marginDegs, spatialSupportSamplesNum, ...
                'forRGCindices', sortedRGCindices(1:min([maxVisualizedRFs numel(sortedRGCindices)])));
        end
    else
        theRetinalRFcenterMaps = [];
    end

    
    if (isempty(xRangeDegs))
        xRangeDegs = obj.sizeDegs(1)*1.05;
    end

    if (isempty(yRangeDegs))
        yRangeDegs = obj.sizeDegs(2)*1.05;
    end

    % Plot part of the input cone mosaic
    if (isempty(xLimsDegs))
        xLimsDegs = mRGCmosaicCenterDegs(1) + 0.5*xRangeDegs*[-1 1];
    end

    if (isempty(yLimsDegs))
        yLimsDegs = mRGCmosaicCenterDegs(2) + 0.5*yRangeDegs*[-1 1];
    end
    
    if (isinf(xLimsDegs)) 
        xLimsDegs = mRGCmosaicCenterDegs(1) + 0.5*obj.inputConeMosaic.sizeDegs(1)*[-1 1];
    end

    if (isinf(yLimsDegs)) 
        yLimsDegs = mRGCmosaicCenterDegs(2) + 0.5*obj.inputConeMosaic.sizeDegs(2)*[-1 1];
    end

    xyLimsDegs = min([xLimsDegs yLimsDegs]);
    if (xyLimsDegs < 0.5)
        xyTicksDegs = 0.05;
    elseif (xyLimsDegs < 1.0)
        xyTicksDegs = 0.2;
    elseif (xyLimsDegs < 2.5)
        xyTicksDegs = 0.2;
    elseif (xyLimsDegs < 5.0)
        xyTicksDegs = 0.5;
    elseif (xyLimsDegs < 10)
        xyTicksDegs = 1.0;
    else
        xyTicksDegs = 2.0;
    end
    
    if (isempty(domainVisualizationTicks))
        xTicks = sign(mRGCmosaicCenterDegs(1)) * round(abs(mRGCmosaicCenterDegs(1)*10))/10 + xyTicksDegs*(-10:1:10);
        yTicks = sign(mRGCmosaicCenterDegs(2)) * round(abs(mRGCmosaicCenterDegs(2)*10))/10 + xyTicksDegs*(-10:1:10);
        domainVisualizationTicks = struct('x', xTicks, 'y', yTicks);
    end

    % Superimpose cone mosaic
    if (superimposedConeMosaic)
        obj.inputConeMosaic.visualize(...
            'figureHandle', figureHandle, ...
            'axesHandle', axesHandle, ...
            'visualizedConeAperture', 'lightCollectingArea4sigma', ...
            'visualizedConeApertureThetaSamples', 20, ...
            'conesAlpha', 0.5, ...
            'withSuperimposedOpticalImage', superimposedOpticalImage, ...
            'withSuperimposedPSF', superimposedPSF, ...
            'domain', 'degrees', ...
            'domainVisualizationLimits', [xLimsDegs(1) xLimsDegs(2) yLimsDegs(1) yLimsDegs(2)], ...
            'domainVisualizationTicks', domainVisualizationTicks, ...
            'clearAxesBeforeDrawing', false, ...
            'noXLabel', noXLabel, ...
            'noYLabel', noYLabel, ...
            'backgroundColor', 'none', ...
            'fontSize', fontSize, ...
            'plotTitle', plotTitle);
   end

  
    if (retinalMeridianAxesLabeling)
        if (~noXLabel)
            if (obj.eccentricityDegs(1) ~= 0)
                % Change the x-label to display the horizontal retinal meridian 
                xlabel(axesHandle, sprintf('%s (degs)', strrep(obj.horizontalRetinalMeridian, 'meridian', 'retina')));
            else
                xlabel(axesHandle, '\leftarrow temporal retina               (degs)                   nasal retina \rightarrow ');
            end
        end

        if (~noYLabel)
            if (obj.eccentricityDegs(2) ~= 0)
                % Change the y-label to display the vertical retinal meridian 
                ylabel(axesHandle, sprintf('%s (degs)', strrep(obj.verticalRetinalMeridian, 'meridian', 'retina')));
            else
                ylabel(axesHandle, '\leftarrow inferior retina                         (degs)                          superior retina\rightarrow  ');
            end
        end
    end

    hold(axesHandle, 'on');

    % Plot the RF maps
    if (isempty(maxVisualizedRFs))
        maxVisualizedRFs = numel(sortedRGCindices);
    end


    if (isempty(inputPoolingVisualization)) && (maxVisualizedRFs == 0)
        plot(axesHandle, obj.rgcRFpositionsDegs(:,1), obj.rgcRFpositionsDegs(:,2), 'k+');
    end

    for iRGC = 1:numel(sortedRGCindices)

        if (iRGC > maxVisualizedRFs)
            continue;
        end

        % Retrieve the RGCindex
        targetRGCindex = sortedRGCindices(iRGC);

        if (~isempty(theRetinalRFcenterMaps))
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
            inputsCentroid = obj.rgcRFpositionsDegs(targetRGCindex,:);
            theCenterConnectedCones = s.inputConeIndices(find(s.inputConeWeights>0.001));
              
            if (numel(s.inputConeIndices)>1) && (showConnectionsToCones)
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
                        'k-', 'LineWidth', contourLineWidth, 'Color', coneColor*0.5);
    
                    plot(axesHandle,[inputsCentroid(1) obj.inputConeMosaic.coneRFpositionsDegs(theCenterConnectedCones(iInput),1)], ...
                        [inputsCentroid(2) obj.inputConeMosaic.coneRFpositionsDegs(theCenterConnectedCones(iInput),2)], ...
                        'k-', 'LineWidth', contourLineWidth*0.5, 'Color', coneColor);
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
                'lineWidth', contourLineWidth, ...
                'edgeAlpha', 0.7);
        else
            switch (inputPoolingVisualization)
                case 'centerOnly'
                    % Retrieve the center cone indices & weights
                    centerConnectivityVector = full(squeeze(obj.rgcRFcenterConePoolingMatrix(:, targetRGCindex)));
                    centerConeIndices = find(centerConnectivityVector > 0.0001);
                    centerConeWeights = reshape(centerConnectivityVector(centerConeIndices), [1 1 numel(centerConeIndices)]);
                    renderConePoolingWeights(axesHandle, ...
                        obj.rgcRFpositionsDegs(targetRGCindex,:), ...
                        obj.inputConeMosaic.coneRFpositionsDegs(centerConeIndices,:), ...
                        centerConeWeights, ...
                        [], [], ...
                        max(centerConeWeights));

                case 'surroundOnly'
                    surroundConnectivityVector = full(squeeze(obj.rgcRFsurroundConePoolingMatrix(:, targetRGCindex)));
                    surroundConeIndices = find(surroundConnectivityVector > 0.0001);
                    surroundConeWeights = reshape(surroundConnectivityVector(surroundConeIndices), [1 1 numel(surroundConeIndices)]);
                    renderConePoolingWeights(axesHandle, ...
                        obj.rgcRFpositionsDegs(targetRGCindex,:), ...
                        [], [], ...
                        obj.inputConeMosaic.coneRFpositionsDegs(surroundConeIndices,:), ...
                        surroundConeWeights, ...
                        max(surroundConeWeights));

                case 'centerAndSurround'
                    centerConnectivityVector = full(squeeze(obj.rgcRFcenterConePoolingMatrix(:, targetRGCindex)));
                    centerConeIndices = find(centerConnectivityVector > 0.0001);
                    centerConeWeights = reshape(centerConnectivityVector(centerConeIndices), [1 1 numel(centerConeIndices)]);
                    surroundConnectivityVector = full(squeeze(obj.rgcRFsurroundConePoolingMatrix(:, targetRGCindex)));
                    surroundConeIndices = find(surroundConnectivityVector > 0.0001);
                    surroundConeWeights = reshape(surroundConnectivityVector(surroundConeIndices), [1 1 numel(surroundConeIndices)]);
                    renderConePoolingWeights(axesHandle, ...
                        obj.rgcRFpositionsDegs(targetRGCindex,:), ...
                        obj.inputConeMosaic.coneRFpositionsDegs(centerConeIndices,:), ...
                        centerConeWeights, ...
                        obj.inputConeMosaic.coneRFpositionsDegs(surroundConeIndices,:), ...
                        surroundConeWeights, ...
                        max(centerConeWeights));
            end
        end
   end % iRGC

   if (~isempty(eccentricitySamplingGrid))
       plot(axesHandle, eccentricitySamplingGrid(:,1), eccentricitySamplingGrid(:,2), 'k+', 'LineWidth', 3.0, 'MarkerSize', 20);
       plot(axesHandle, eccentricitySamplingGrid(:,1), eccentricitySamplingGrid(:,2), 'c+', 'LineWidth', 1.0, 'MarkerSize', 16);
   end

   set(axesHandle, 'XLim', [xLimsDegs(1) xLimsDegs(2)], 'YLim', [yLimsDegs(1) yLimsDegs(2)]);
   drawnow;
end

function renderConePoolingWeights(ax, rgcRFpositionsDegs, rfCenterConePositionsDegs, rfCenterConeWeights, ...
                  rfSurroundConePositionsDegs, rfSurroundConeWeights, maxConeWeight)

    for iCone = 1:numel(rfSurroundConeWeights)
        if (~isempty(rfCenterConeWeights))
            surroundBoostFactor = 20;
        else
            surroundBoostFactor = 1;
        end
        w = 0.5 + surroundBoostFactor * rfSurroundConeWeights(iCone)/maxConeWeight*7;
        plot(ax, [rgcRFpositionsDegs(1), rfSurroundConePositionsDegs(iCone,1)], ...
                 [rgcRFpositionsDegs(2), rfSurroundConePositionsDegs(iCone,2)], ...
                 'k-', 'LineWidth', w);
    end
    
    if (numel(rfCenterConeWeights) == 1)
        plot(ax, rfCenterConePositionsDegs(1,1), rfCenterConePositionsDegs(1,2), ...
            'wo', 'LineWidth', 2.0, 'MarkerSize', 20);
    else
        for iCone = 1:numel(rfCenterConeWeights)
            w = 0.5 + rfCenterConeWeights(iCone)/maxConeWeight*7;
            plot(ax, [rgcRFpositionsDegs(1), rfCenterConePositionsDegs(iCone,1)], ...
                     [rgcRFpositionsDegs(2), rfCenterConePositionsDegs(iCone,2)], ...
                     'w-', 'LineWidth', w);
        end
    end
    

end

