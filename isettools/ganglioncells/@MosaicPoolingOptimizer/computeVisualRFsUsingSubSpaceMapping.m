function computeVisualRFsUsingSubSpaceMapping(...
            theComputeReadyMRGCmosaic, opticsToEmploy, ...
            stimSizeDegs, stimPositionDegs, ...
            maxSFLimit, maxSFToBeAnalyzed, rfMappingPixelMagnificationFactor, ...
            stimulusChromaticity, coneFundamentalsOptimizedForStimPosition, ...
            coneMosaicResponsesFileName, ...
            mRGCMosaicResponsesFileName, ...
            optimallyMappedRFmapsFileName, ...
            reComputeInputConeMosaicResponses, ...
            reComputeMRGCMosaicResponses, ...
            reComputeRFs, ...
            visualizeOptimallyMappedRFmapLocations, ...
            varargin)

    p = inputParser;
    p.addParameter('parPoolSize', [], @(x)(isempty(x)||(isscalar(x))));
    p.addParameter('visualizedResponses', false, @islogical);
    p.addParameter('msequencePixelSizeDegs', 0.01, @isscalar);
    p.addParameter('visualizedRGCindex', [], @(x)(isempty(x)||(isscalar(x))));
    p.parse(varargin{:});
    parPoolSize = p.Results.parPoolSize;
    visualizedResponses = p.Results.visualizedResponses;
    msequencePixelSizeDegs = p.Results.msequencePixelSizeDegs;
    visualizedRGCindex = p.Results.visualizedRGCindex;

    if (isempty(stimPositionDegs))
        stimPositionDegs = theComputeReadyMRGCmosaic.eccentricityDegs;
    end

    % Encode examined spatial position, pixel magnification and maxSF
    if (isempty(maxSFLimit))
        positionPostFix = sprintf('_atPosition_%2.2f_%2.2f_PixMag_%2.1f_maxSF_OPTIMAL.mat', stimPositionDegs(1), stimPositionDegs(2), rfMappingPixelMagnificationFactor);
    else
        positionPostFix = sprintf('_atPosition_%2.2f_%2.2f_PixMag_%2.1f_maxSF_%2.0fCPD.mat', stimPositionDegs(1), stimPositionDegs(2), rfMappingPixelMagnificationFactor, maxSFLimit);
    end

    coneMosaicResponsesFileName = strrep(coneMosaicResponsesFileName, '.mat', positionPostFix);
    mRGCMosaicResponsesFileName = strrep(mRGCMosaicResponsesFileName, '.mat', positionPostFix);
    optimallyMappedRFmapsFileName = strrep(optimallyMappedRFmapsFileName, '.mat', positionPostFix);

    if (reComputeInputConeMosaicResponses || reComputeMRGCMosaicResponses || reComputeRFs)

        % Compute the RF maps for ALL cells using stimuli at this 
        % position. Only some cells will be optimally mapped at each
        % position. These optimally derived RF maps are extracted by
        % the optimalyMappedRFsAtMosaicPosition() function
        computeRFmapsForAllCellsUsingStimuliAtTargetPosition(...
            theComputeReadyMRGCmosaic, opticsToEmploy, ...
            stimSizeDegs, stimPositionDegs, ...
            maxSFLimit, maxSFToBeAnalyzed, rfMappingPixelMagnificationFactor, ...
            stimulusChromaticity, coneFundamentalsOptimizedForStimPosition, ...     
            parPoolSize, ...
            coneMosaicResponsesFileName, ...
            mRGCMosaicResponsesFileName, ...
            optimallyMappedRFmapsFileName, ...
            reComputeInputConeMosaicResponses, ...
            reComputeMRGCMosaicResponses, ...
            reComputeRFs, ...
            visualizedResponses);
    end

    if (visualizeOptimallyMappedRFmapLocations)
        visualizeAllOptimallyMappedRFmapLocations(optimallyMappedRFmapsFileName, ...
            theComputeReadyMRGCmosaic, visualizedRGCindex, stimulusChromaticity, msequencePixelSizeDegs);
    end

end

function visualizeAllOptimallyMappedRFmapLocations(optimallyMappedRFmapsFileName,...
    theMRGCMosaic, visualizedRGCindex, stimulusChromaticity, msequencePixelSizeDegs)

    % Load all the optimally mapped visual RF maps
    fprintf('Loading optimally mapped subspace RF maps from %s\n', optimallyMappedRFmapsFileName);
    load(optimallyMappedRFmapsFileName, 'optimallyMappedVisualRFmaps', 'indicesOfOptimallyMappedRGCs');

    hFig = figure(22); clf;
    set(hFig, 'Color', [1 1 1], 'Position', [10 10 1200 1280]);

    subplotPosVectors = NicePlot.getSubPlotPosVectors(...
       'rowsNum', 2, ...
       'colsNum', 2, ...
       'heightMargin',  0.06, ...
       'widthMargin',    0.05, ...
       'leftMargin',     0.06, ...
       'rightMargin',    0.05, ...
       'bottomMargin',   0.05, ...
       'topMargin',      0.02);


    for iCell = 1:numel(indicesOfOptimallyMappedRGCs)
        theRGCindex = indicesOfOptimallyMappedRGCs(iCell);

        if ((~isempty(visualizedRGCindex)) && (theRGCindex ~= visualizedRGCindex))
            continue;
        end

        [~,~,~, theCenterConeTypes, theCenterConeIndices] = ...
            theMRGCMosaic.centerConeTypeWeights(theRGCindex);

        % The RF center position on the retina
        theRGCRetinalRFcenterPositionDegs = theMRGCMosaic.rgcRFpositionsDegs(theRGCindex,:);

        [~, ~, theSurroundConeTypes, theSurroundConeIndices] = ...
            theMRGCMosaic.surroundConeTypeWeights(theRGCindex, theCenterConeIndices);

        idx = find(theSurroundConeTypes == cMosaic.LCONE_ID);
        surroundLconeIndices = theSurroundConeIndices(idx);

        idx = find(theSurroundConeTypes == cMosaic.MCONE_ID);
        surroundMconeIndices = theSurroundConeIndices(idx);

        d = optimallyMappedVisualRFmaps{iCell};
        d = generateMsequenceRFmap(d, msequencePixelSizeDegs);

        xLim = [d.spatialSupportDegsX(1) d.spatialSupportDegsX(end)];
        %xLim = theRGCRetinalRFcenterPositionDegs(1) + 0.08*[-1 1];

        yLim = [d.spatialSupportDegsY(1) d.spatialSupportDegsY(end)];
        %yLim = theRGCRetinalRFcenterPositionDegs(2) + 0.08*[-1 1];

        % The RF center position in the visual field
        [~,idx] = max(abs(d.theRFmap(:)));
        [row,col] = ind2sub(size(d.theRFmap), idx);
        theRGCVisualRFcenterPositionDegs = [d.spatialSupportDegsX(col) d.spatialSupportDegsY(row)];

        ax = subplot('Position', subplotPosVectors(1,1).v);
        renderSubspaceRFmap(ax, d, ...
            theMRGCMosaic.inputConeMosaic.coneRFpositionsDegs(theCenterConeIndices,:), ...
            theMRGCMosaic.inputConeMosaic.coneRFpositionsDegs(surroundLconeIndices,:), ...
            theMRGCMosaic.inputConeMosaic.coneRFpositionsDegs(surroundMconeIndices,:), ...
            theRGCindex, theRGCRetinalRFcenterPositionDegs, theRGCVisualRFcenterPositionDegs, ...
            stimulusChromaticity, xLim, yLim);

        xProfile = sum(d.theRFmap, 1);
        yProfile = sum(d.theRFmap, 2);
        maxP = max([max(abs(xProfile(:))) max(abs(yProfile(:)))]);

        ax = subplot('Position', subplotPosVectors(2,2).v);
        renderXProfile(ax, d.spatialSupportDegsX, xProfile, maxP, theRGCVisualRFcenterPositionDegs, xLim);

        ax = subplot('Position', subplotPosVectors(1,2).v);
        %renderYProfile(ax, d.spatialSupportDegsY, yProfile, maxP, theRGCVisualRFcenterPositionDegs, yLim);
        renderSFTuningPlot(ax, d.spatialFrequencySupport, d.theSFtuningMap);

        ax = subplot('Position', subplotPosVectors(2,1).v);
        renderMsequenceRFmap(ax, d, theRGCVisualRFcenterPositionDegs, msequencePixelSizeDegs, xLim, yLim);

        drawnow;
    end

end

function d = generateMsequenceRFmap(d, msequencePixelSizeDegs) 
    dx = d.spatialSupportDegsX(2)-d.spatialSupportDegsX(1);
    msequencePixelSize = round(msequencePixelSizeDegs / dx);
    msequencePixelSizeDegs = dx*msequencePixelSize;

    msequencePixelKernel = ones(msequencePixelSize,msequencePixelSize);
    msequencePixelsNum = round((d.spatialSupportDegsX(end)-d.spatialSupportDegsX(1))/msequencePixelSizeDegs);

    theMsequenceMap = conv2(d.theRFmap, msequencePixelKernel, 'same');

    [X,Y] = meshgrid(d.spatialSupportDegsX, d.spatialSupportDegsY);
    F = scatteredInterpolant(double(X(:)), double(Y(:)), double(theMsequenceMap(:)));

    msequenceSpatialSupport = (1:msequencePixelsNum)*msequencePixelSizeDegs;
    msequenceSpatialSupport = msequenceSpatialSupport - mean(msequenceSpatialSupport);
    d.theMsequenceRFmapSpatialSupportX = mean(d.spatialSupportDegsX) + msequenceSpatialSupport;
    d.theMsequenceRFmapSpatialSupportY = mean(d.spatialSupportDegsY) + msequenceSpatialSupport;

    [X,Y] = meshgrid(d.theMsequenceRFmapSpatialSupportX, d.theMsequenceRFmapSpatialSupportY);
    d.theMsequenceRFmap = reshape(F(X(:),Y(:)), [msequencePixelsNum msequencePixelsNum]);

end


function shadedAreaPlotX(ax,x,y, baseline, faceColor, edgeColor, faceAlpha, lineWidth)
    x = [x(:); flipud(x(:))];
    y = [y(:); y(:)*0+baseline];
    px = reshape(x, [1 numel(x)]);
    py = reshape(y, [1 numel(y)]);
    pz = -10*eps*ones(size(py)); 
    patch(ax,px,py,pz,'FaceColor',faceColor,'EdgeColor', edgeColor, 'FaceAlpha', faceAlpha, 'LineWidth', lineWidth);
end

function shadedAreaPlotY(ax, x,y, baseline, faceColor, edgeColor, faceAlpha, lineWidth)
    xx = [y(:); flipud(y(:))*0+baseline];
    yy = [x(:); flipud(x(:))];
    px = reshape(xx, [1 numel(xx)]);
    py = reshape(yy, [1 numel(yy)]);
    pz = -10*eps*ones(size(py)); 
    patch(ax,px,py,pz,'FaceColor',faceColor,'EdgeColor', edgeColor, 'FaceAlpha', faceAlpha, 'LineWidth', lineWidth);
end

function renderXProfile(ax, spatialSupport, profile, maxP, theRGCRFcenterPositionDegs, xLim)

        cla(ax)

        shadedAreaPlotX(ax, spatialSupport, profile, 0, [0.5 0.5 0.5], [0.5 0.5 0.5], 0.5, 1.0);
        hold(ax, 'on');
        plot(ax, spatialSupport, profile*0, 'k-');

        plot(theRGCRFcenterPositionDegs(1)*[1 1], maxP*[-1 1], 'b-', 'LineWidth', 1.0);
        hold(ax, 'off');

        axis(ax, 'square')
        set(ax, 'XLim', xLim);
        set(ax, 'YLim', maxP*[-1 1]);
        set(ax, 'XTick', -10:0.05:10);
        set(ax, 'XColor', [0 0 0], 'YColor', [1 0 0]);
        xlabel(ax, 'space, x (degs)');
        set(ax, 'FontSize', 16)
end

function renderYProfile(ax, spatialSupport, profile, maxP, theRGCRFcenterPositionDegs, yLim)
    cla(ax)
    yyaxis(ax, 'left');
    set(ax, 'YColor', 'none');
    yyaxis(ax, 'right');

    shadedAreaPlotY(ax, spatialSupport, profile, 0, [0.5 0.5 0.5], [0.5 0.5 0.5], 0.5, 1.0);
    hold(ax, 'on');
    plot(ax, profile*0, spatialSupport, 'k-');

    plot(maxP*[-1 1], theRGCRFcenterPositionDegs(2)*[1 1], 'b-', 'LineWidth', 1.0);
    hold(ax, 'off');

    axis(ax, 'square');
    set(ax, 'YLim', yLim);
    set(ax, 'YTick', -10:0.05:10);
    set(ax, 'XDir', 'reverse', 'XLim', maxP*[-1 1]);
    set(ax, 'YColor', [0 0 0], 'XColor', [1 0 0]);
    ylabel(ax, 'space, y (degs)');
    set(ax, 'FontSize', 16);
end


function renderSubspaceRFmap(ax, d, centerConePositions, surroundLconePositions, surroundMconePositions, theRGCindex, ...
    theRGCRFcenterPositionDegs, theVisualRGCRFcenterPositionDegs, stimulusChromaticity, xLim, yLim)

    imagesc(ax, d.spatialSupportDegsX,  d.spatialSupportDegsY, d.theRFmap);
    hold(ax, 'on');    

    
    scatter(ax, surroundLconePositions(:,1), surroundLconePositions(:,2), ...
         16*16, 'MarkerFaceColor', [1 0.5 0.5], 'MarkerEdgeColor', [1 0 0], 'MarkerFaceAlpha', 0.5, 'LineWidth', 1.0);
        
    scatter(ax, surroundMconePositions(:,1), surroundMconePositions(:,2), ...
         16*16, 'MarkerFaceColor', [0.5 0.8 0.5], 'MarkerEdgeColor', [0 0.8 0], 'MarkerFaceAlpha', 0.5, 'LineWidth', 1.0);
    
    plot(ax, centerConePositions(:,1), centerConePositions(:,2), 'ko', 'MarkerSize', 16, 'LineWidth', 2.0);

    %plot(ax, theRGCRFcenterPositionDegs(1), theRGCRFcenterPositionDegs(2), 'yx', 'MarkerSize', 24, 'LineWidth', 2.0);
    plot(theVisualRGCRFcenterPositionDegs(1)*[1 1], [-100 100], 'b-', 'LineWidth', 1.0);
    plot([-100 100], theVisualRGCRFcenterPositionDegs(2)*[1 1], 'b-', 'LineWidth', 1.0);

    hold(ax, 'off')
    axis(ax,'image'); axis(ax,'xy');
    set(ax, 'CLim', max(abs(d.theRFmap(:)))*[-1 1], ...
            'XLim', xLim, ...
            'YLim', yLim, ...
            'FontSize', 16 ...
    );
    set(ax, 'XTick', -10:0.05:10);
    set(ax, 'YTick', -10:0.05:10);
    cLUT = brewermap(1024, '*RdBu');
    colormap(ax, cLUT);
    set(ax, 'Color', cLUT(512,:));
    xlabel(ax, 'space, x (degs)');
    ylabel(ax, 'space, y (degs)');
    title(ax, sprintf('RGC %d - %s RF ([%2.2f ... %2.2f])', theRGCindex, stimulusChromaticity, min(d.theRFmap(:)), max(d.theRFmap(:))));
end


function renderMsequenceRFmap(ax, d, theRGCRFcenterPositionDegs, msequencePixelSizeDegs, xLim, yLim)
    imagesc(ax, d.theMsequenceRFmapSpatialSupportX, d.theMsequenceRFmapSpatialSupportY, d.theMsequenceRFmap);
    hold(ax, 'on');
    zLevelsNeg = 0.05: 0.1: 1;
    zLevelsPos = 0.05: 0.1: 1;
    contour(ax, d.theMsequenceRFmapSpatialSupportX,  d.theMsequenceRFmapSpatialSupportY, d.theMsequenceRFmap, max(abs(d.theMsequenceRFmap(:)))*zLevelsPos, ...
        'LineStyle', '-', 'LineWidth', 1.0, 'EdgeColor', [0 0 0]);
    contour(ax, d.theMsequenceRFmapSpatialSupportX,  d.theMsequenceRFmapSpatialSupportY, d.theMsequenceRFmap, -max(abs(d.theMsequenceRFmap(:)))*zLevelsNeg, ...
        'LineStyle', '--', 'LineWidth', 1.0, 'EdgeColor', [0 0 0]);
    plot(theRGCRFcenterPositionDegs(1)*[1 1], [-100 100], 'b-', 'LineWidth', 1.0);
    plot([-100 100], theRGCRFcenterPositionDegs(2)*[1 1], 'b-', 'LineWidth', 1.0);
    hold(ax, 'off');
    axis(ax,'image'); axis(ax,'xy');
    set(ax, 'CLim', max(abs(d.theMsequenceRFmap(:)))*[-1 1], ...
            'XLim', xLim, ...
            'YLim', yLim, ...
            'FontSize', 16 ...
    );
    set(ax, 'XTick', -10:0.05:10);
    set(ax, 'YTick', -10:0.05:10);
    cLUT = brewermap(1024, '*RdBu');
    colormap(ax, cLUT);
    set(ax, 'Color', cLUT(512,:));
    xlabel(ax, 'space, x (degs)');
    ylabel(ax, 'space, y (degs)');
    %title(ax, sprintf('m-sequence mapped RF (pixel size: %2.3f arc min)', msequencePixelSizeDegs*60));
    title(ax, sprintf('m-sequence mapped RF'));
end


function renderSFTuningPlot(ax, spatialFrequencySupport, theSFtuningMap)
    
    % Take the magnitude
    omega = (numel(spatialFrequencySupport)-1)/2;
    
    theSFtuningMapMag = 0*theSFtuningMap;
    for iy = -omega:omega
        for ix = -omega:omega
            ixNeg = -ix; iyNeg = -iy;
            theSFtuningMapMag(iy+1+omega,ix+1+omega) = ...
                sqrt(theSFtuningMap(iy+1+omega,ix+1+omega).^2 + theSFtuningMap(iyNeg+1+omega, ixNeg+1+omega).^2);
        end
    end

    linearScale = false;
    logScale = false;

    radialSlices = true;

    if (radialSlices)
        deltaAngle = 15;
        idx = find(spatialFrequencySupport >=0);
        midPoint = find(spatialFrequencySupport == 0);
        hiResSFaxis = linspace(0,max(spatialFrequencySupport(:)), 100);
        anglesExamined = 0:deltaAngle:(360-deltaAngle);
        lineColors = brewermap(numel(anglesExamined), 'Set1');
        theLegends = {};
        for iAngle = 1:numel(anglesExamined)
            theSlice = imrotate(theSFtuningMapMag, anglesExamined(iAngle),'bilinear','crop');
            plot(ax, 0.1+hiResSFaxis, interp1(spatialFrequencySupport(idx), theSlice(midPoint,idx), hiResSFaxis), ...
                '-', 'LineWidth', 2.0, 'Color', [0 0 0]);
            hold(ax, 'on');
            plot(ax, 0.1+hiResSFaxis, interp1(spatialFrequencySupport(idx), theSlice(midPoint,idx), hiResSFaxis), ...
                '-', 'LineWidth', 1.5, 'Color', squeeze(lineColors(iAngle,:)));
            theLegends{numel(theLegends)+1} = sprintf('%d^o', anglesExamined(iAngle));
        end
        legend(ax, theLegends, 'NumColumns', 2, 'Location', 'SouthWest');
        hold(ax, 'off')
        axis(ax, 'square');
        set(ax, 'XScale', 'log', 'XLim', [0.5 100], ...
                    'YLim', [0 max(theSFtuningMapMag(:))], ...
                    'YTick', 0:0.02:1, ...
                'FontSize', 16 ...
            );
        grid(ax, 'on'); box(ax, 'off');
            set(ax, 'XTick', [0.1 0.3 1 3 10 30 100], 'YTickLabel', {});
            xlabel(ax, 'spatial frequency, x (cpd)');
            ylabel(ax, 'STF');
        title(ax, 'STF (slices through 2DFT (RFmap)')
    end

    if (linearScale)
        imagesc(ax, spatialFrequencySupport, spatialFrequencySupport, theSFtuningMapMag);
        hold(ax, 'on')
        plot([0 0], [-100 100], 'c-', 'LineWidth', 1.0);
        plot([-100 100], [0 0],'c-', 'LineWidth', 1.0);
        hold(ax, 'off');
        axis(ax,'image'); axis(ax,'xy');
        set(ax, 'XLim', [-90 90], ...
                'YLim', [-90 90], ...
                'FontSize', 16 ...
        );
        set(ax, 'XTick', -100:10:100);
        set(ax, 'YTick', -100:10:100);
        set(ax, 'CLim', [0 max(theSFtuningMapMag(:))]);
        colormap(ax,gray(1024));
        xlabel(ax, 'spatial frequency, x (cpd)');
        ylabel(ax, 'spatial frequency, y (cpd)');
    end

    if (logScale)
        nPts = 20;
        minSF = spatialFrequencySupport(2)-spatialFrequencySupport(1);
        sfLogScalingPos = logspace(log10(minSF), log10(100), nPts);
        sfLogScalingFull = [fliplr(-sfLogScalingPos) 0 sfLogScalingPos];
        
        [X,Y] = ndgrid(spatialFrequencySupport, spatialFrequencySupport);
        F = griddedInterpolant(X,Y,theSFtuningMapMag);
        [X,Y] = ndgrid(sfLogScalingFull,sfLogScalingFull);
        theSFtuningMapMagLogScale = reshape(F(X, Y), [numel(sfLogScalingFull) numel(sfLogScalingFull)]);


        imagesc(ax, -nPts:nPts, -nPts:nPts, theSFtuningMapMagLogScale);
        hold(ax, 'on')
        plot([0 0], [-100 100], 'c-', 'LineWidth', 1.0);
        plot([-100 100], [0 0],'c-', 'LineWidth', 1.0);
        hold(ax, 'off');
        axis(ax,'image'); axis(ax,'xy');
        set(ax, 'XLim', [ -nPts  nPts], ...
                'YLim', [ -nPts  nPts], ...
                'FontSize', 16 ...
        );
        set(ax, 'XTick', -nPts:2:nPts, 'XTickLabel', sprintf('%2.0f\n', sfLogScalingFull(1:2:end)));
        set(ax, 'YTick', -nPts:2:nPts, 'YTickLabel', sprintf('%2.0f\n', sfLogScalingFull(1:2:end)));
        set(ax, 'CLim', [0 max(theSFtuningMapMagLogScale(:))]);
        xtickangle(ax, 0);
        colormap(ax,gray(1024));
        xlabel(ax, 'spatial frequency, x (cpd)');
        ylabel(ax, 'spatial frequency, y (cpd)');
        
    end


end


function computeRFmapsForAllCellsUsingStimuliAtTargetPosition( ...
            theComputeReadyMRGCmosaic, opticsToEmploy, ...
            stimSizeDegs, stimPositionDegs, ...
            maxSFLimit, maxSFToBeAnalyzed, rfMappingPixelMagnificationFactor, ...
            stimulusChromaticity, coneFundamentalsOptimizedForStimPosition, ...
            parPoolSize,...
            coneMosaicResponsesFileName, ...
            mRGCMosaicResponsesFileName, ...
            optimallyMappedRFmapsFileName, ...
            reComputeInputConeMosaicResponses, ...
            reComputeMRGCMosaicResponses, ...
            reComputeRFs, ...
            visualizedResponses)

    % Compute responses of the input cone mosaic to the subspace RF mapping stimuli
    if (reComputeInputConeMosaicResponses)
        fprintf('Cone mosaic subspace responses will be saved to %s \n', coneMosaicResponsesFileName);
        MosaicPoolingOptimizer.generateInputConeMosaicSubspaceRFmappingLinearResponses(...
            theComputeReadyMRGCmosaic, opticsToEmploy, ...
            stimSizeDegs, stimPositionDegs, ...
            stimulusChromaticity, coneFundamentalsOptimizedForStimPosition, ...
            coneMosaicResponsesFileName, ...
            'maxSFLimit', maxSFLimit, ...
            'rfMappingPixelMagnificationFactor', rfMappingPixelMagnificationFactor, ...
            'visualizedResponses', visualizedResponses, ...
            'parPoolSize', parPoolSize);
    end


    if (reComputeMRGCMosaicResponses)
        fprintf('\nLoading cone mosaic subspace modulation responses and Hartley spatial modulation patterns ...');
        % Load the previously computed responses
        load(coneMosaicResponsesFileName, ...
            'HartleySpatialModulationPatterns', 'spatialSupportDegs', 'stimParams', 'lIndices', 'mIndices', ...
            'theConeMosaicSubspaceLinearModulationResponses');
        HartleySpatialModulationPatterns = single(HartleySpatialModulationPatterns);
        fprintf('Done loading !\n');

        [HartleyStimNum, nCones] = size(theConeMosaicSubspaceLinearModulationResponses);

        fprintf('MRGC mosaic subspace RF maps and responses will be saved to %s \n', mRGCMosaicResponsesFileName);
        fprintf('Computing visual subspace RF mapping responses for all RGCs in the mosaic ... \n');
        nTimePoints = 1;
        nTrials = 1;
        theConeMosaicResponseTemporalSupportSeconds = [0];
        
        theMRGCMosaicSubspaceRFmappingLinearResponses = zeros(...
                     HartleyStimNum, theComputeReadyMRGCmosaic.rgcsNum, 'single');

        % Use all processors
        [shutdownParPoolOnceCompleted, numWorkers] = MosaicPoolingOptimizer.resetParPool([]);
    
        % Noise-free responses
        theComputeReadyMRGCmosaic.noiseFlag = 'none';

        parfor iStim = 1:size(theConeMosaicSubspaceLinearModulationResponses,1)
             fprintf('Computing mRGC mosaic response for Hartley pattern %d of %d (using %d parallel processes).\n', ...
                 iStim, HartleyStimNum, numWorkers);
             theConeMosaicModulationResponse = squeeze(theConeMosaicSubspaceLinearModulationResponses(iStim,:));
             theConeMosaicModulationResponse = reshape(theConeMosaicModulationResponse, [nTrials nTimePoints nCones]);

             % Compute mRGC mosaic responses based on the cone mosaic modulation responses
             [theMRGCMosaicResponse, ~, theMRGCresponseTemporalSupportSeconds] = ...
                    theComputeReadyMRGCmosaic.compute(theConeMosaicModulationResponse, theConeMosaicResponseTemporalSupportSeconds);
             theMRGCMosaicSubspaceRFmappingLinearResponses(iStim,:) = single(squeeze(theMRGCMosaicResponse(nTrials, nTimePoints ,:)));
        end
    
        fprintf('\nSaving computed mRGCRF mosaic SUBSPACE RF mapping linear responses to %s ...', mRGCMosaicResponsesFileName);
        save(mRGCMosaicResponsesFileName, ...
            'stimParams', ...
            'theMRGCMosaicSubspaceRFmappingLinearResponses', ...
            'spatialSupportDegs', 'lIndices', 'mIndices', ...
            '-v7.3');
    end

    if (reComputeRFs)
        load(coneMosaicResponsesFileName, 'HartleySpatialModulationPatterns', 'spatialSupportDegs', 'stimParams', 'lIndices', 'mIndices');
        HartleySpatialModulationPatterns = single(HartleySpatialModulationPatterns);

        % Load theMRGCMosaicSubspaceRFmappingLinearResponses
        load(mRGCMosaicResponsesFileName, ...
                'theMRGCMosaicSubspaceRFmappingLinearResponses');
    
        % Determine indices of RGCs whose RF lie within the stimulus region
        indicesOfOptimallyMappedRGCs = MosaicPoolingOptimizer.indicesOfOptimallyMappedRGCsAtThisPosition(theComputeReadyMRGCmosaic, ...
            stimPositionDegs, stimSizeDegs);
        

        % Compute RF maps of all cells within the stimulus region
        [theMRGCMosaicOptimallyMappedVisualRFmaps, theMRGCMosaicOptimallyMappedFrequencyTuningMaps, spatialFrequencySupport] = computeRFs(...
            indicesOfOptimallyMappedRGCs, ...
            theMRGCMosaicSubspaceRFmappingLinearResponses, ...
            HartleySpatialModulationPatterns, ...
            spatialSupportDegs, lIndices, mIndices, maxSFToBeAnalyzed);
    
        fprintf('\nSaving computed visual RFs to %s ...', mRGCMosaicResponsesFileName);
        save(mRGCMosaicResponsesFileName, ...
            'theMRGCMosaicOptimallyMappedVisualRFmaps', ...
            'theMRGCMosaicOptimallyMappedFrequencyTuningMaps', ...
            'spatialFrequencySupport', ...
            'indicesOfOptimallyMappedRGCs', '-append');
        fprintf('Done saving! \n');

        % Export visual RF maps for cells that are optimally mapped at this  position
        exportOptimalyMappedRFmaps(stimPositionDegs, mRGCMosaicResponsesFileName, optimallyMappedRFmapsFileName);
    end

end


function [theRFmaps, theFrequencyTuningMaps, spatialFrequencySupport] = computeRFs( ...
    indicesOfOptimallyMappedRGCs, ...
    theSubspaceRFmappingLinearResponses, ...
    HartleySpatialModulationPatterns, ...
    spatialSupportDegs, lIndices, mIndices, maxSFtoBeAnalyzed)

    nStim = size(theSubspaceRFmappingLinearResponses,1);
    cellsNum = size(theSubspaceRFmappingLinearResponses,2);
    pixelsNum = size(HartleySpatialModulationPatterns,2);

    m = max(abs(theSubspaceRFmappingLinearResponses),[],1);
    cellsWithNonZeroResponse = find(m > 0);
    theRFmaps = cell(cellsNum, 1);
    theFrequencyTuningMaps = cell(cellsNum, 1);

    visualizeStimulusSet = true;
    [~,~,~,~, spatialFrequencySupport] = rfMappingStimulusGenerator.HartleySFmap(HartleySpatialModulationPatterns, ...
        spatialSupportDegs, lIndices, mIndices, 1, visualizeStimulusSet, maxSFtoBeAnalyzed);

    parfor iCell = 1:numel(indicesOfOptimallyMappedRGCs)
        fprintf('Computing visual RF by accumulating Hartley patterns for the %d of %d optimally mapped RGC... \n', iCell, numel(indicesOfOptimallyMappedRGCs));
        theRGCindex = indicesOfOptimallyMappedRGCs(iCell);

        if (ismember(theRGCindex, cellsWithNonZeroResponse))
            theRFmap = zeros(pixelsNum, pixelsNum, 'single');
            theFrequencyTuningMap = zeros(numel(spatialFrequencySupport), numel(spatialFrequencySupport));
            allResponses = squeeze(theSubspaceRFmappingLinearResponses(:,theRGCindex));

            for iStim = 1:nStim
                [theStimulusFxCoord, theStimulusFyCoord, stimulusInInsideAnalyzedSFregion, theStimulusFrequency] = rfMappingStimulusGenerator.HartleySFmap(...
                    HartleySpatialModulationPatterns, spatialSupportDegs, lIndices, mIndices, iStim, false, maxSFtoBeAnalyzed);

                if (stimulusInInsideAnalyzedSFregion)
                    theFrequencyTuningMap(theStimulusFyCoord, theStimulusFxCoord) = allResponses(iStim);
                    theRFmap = theRFmap + ...
                            single(squeeze(HartleySpatialModulationPatterns(iStim,:,:)) * allResponses(iStim));
                end
            end
            theRFmaps{iCell} = theRFmap;
            theFrequencyTuningMaps{iCell} = theFrequencyTuningMap;
        end
    end
end


function exportOptimalyMappedRFmaps(stimPositionDegs, mRGCMosaicResponsesFileName, optimallyMappedRFmapsFileName)

    % Check to see if the subspace responses have been generated at this
    % grid position
    fileExists = isfile(mRGCMosaicResponsesFileName);
    if (~fileExists)
        error('File %s does not exist.\n Not extracting optimally mapped visual RF maps at (x,y)=(%2.1f,%2.1f)\n', ...
            mRGCMosaicResponsesFileName, ...
            stimPositionDegs(1), stimPositionDegs(2));
    end

    % Check to see if subspace RF maps have been computed at this grid position
    s = whos('-file', mRGCMosaicResponsesFileName);
    computedVariableNames = cell(1, numel(s));
    for i = 1:numel(s)
        computedVariableNames{i} = s(i).name;
    end
    
    if (~ismember('theMRGCMosaicOptimallyMappedVisualRFmaps', computedVariableNames))
        error('Subspace responses are computed, however RF maps have not been computed yet.\nNot extracting optimally mapped visual RF maps at (x,y) = (%2.1f,%2.1f)\n', ...
            stimPositionDegs(1), stimPositionDegs(2));
    end
    
    fprintf('Loading RF map data from %s. Please wait ...\n', mRGCMosaicResponsesFileName);
    % All good. Extract the computed subspace RF maps
    load(mRGCMosaicResponsesFileName, ...
        'spatialSupportDegs', ...
        'indicesOfOptimallyMappedRGCs', ...
        'theMRGCMosaicOptimallyMappedVisualRFmaps', ...
        'theMRGCMosaicOptimallyMappedFrequencyTuningMaps', ...
        'spatialFrequencySupport');


    optimallyMappedVisualRFmaps = cell(1, numel(indicesOfOptimallyMappedRGCs));

    dx = spatialSupportDegs(2)-spatialSupportDegs(1);
    for iCell = 1:numel(theMRGCMosaicOptimallyMappedVisualRFmaps)
        optimallyMappedVisualRFmaps{iCell} = struct(...
            'theRFmap', theMRGCMosaicOptimallyMappedVisualRFmaps{iCell}, ...
            'spatialSupportDegsX', spatialSupportDegs+stimPositionDegs(1)-dx, ...
            'spatialSupportDegsY', spatialSupportDegs+stimPositionDegs(2)-dx, ...
            'theSFtuningMap', theMRGCMosaicOptimallyMappedFrequencyTuningMaps{iCell}, ...
            'spatialFrequencySupport', spatialFrequencySupport ...
            );
    end


    % Save all the optimally mapped visual RF maps
    fprintf('Saving optimally mapped subspace RF maps to %s\n', optimallyMappedRFmapsFileName);
    save(optimallyMappedRFmapsFileName, ...
        'optimallyMappedVisualRFmaps', 'indicesOfOptimallyMappedRGCs', '-v7.3');

end