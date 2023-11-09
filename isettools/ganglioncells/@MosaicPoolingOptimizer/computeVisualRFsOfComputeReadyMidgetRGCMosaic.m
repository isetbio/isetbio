function computeVisualRFsOfComputeReadyMidgetRGCMosaic(...
            theComputeReadyMRGCmosaic, opticsToEmploy, ...
            stimSizeDegs, stimPositionDegs, stimulusChromaticity,...
            employConeFundamentalsDerivedFromInputConeMosaicAtStimPosition, ...
            maxSFLimit, rfMappingPixelMagnificationFactor, ...
            coneMosaicSubspaceResponsesFileName, ...
            mRGCMosaicSubspaceResponsesFileName, ...
            optimallyMappedSubspaceRFmapsFileName, ...
            reComputeInputConeMosaicSubspaceRFmappingResponses, ...
            reComputeMRGCMosaicSubspaceRFmappingResponses, ...
            reComputeRFs, ...
            visualizeOptimallyMappedRFmapLocations, ...
            varargin)

    p = inputParser;
    p.addParameter('parPoolSize', [], @(x)(isempty(x)||(isscalar(x))));
    p.addParameter('visualizedResponses', false, @islogical);
    p.addParameter('msequencePixelSizeDegs', 0.01, @isscalar);
    p.parse(varargin{:});
    parPoolSize = p.Results.parPoolSize;
    visualizedResponses = p.Results.visualizedResponses;
    msequencePixelSizeDegs = p.Results.msequencePixelSizeDegs;

    if (isempty(stimPositionDegs))
        stimPositionDegs = theComputeReadyMRGCmosaic.eccentricityDegs;
    end

    % Encode examined spatial position, pixel magnification and maxSF
    if (isempty(maxSFLimit))
        positionPostFix = sprintf('_atPosition_%2.2f_%2.2f_PixelMagnification_%2.3f_maxSF_OPTIMAL.mat', stimPositionDegs(1), stimPositionDegs(2), rfMappingPixelMagnificationFactor);
    else
        positionPostFix = sprintf('_atPosition_%2.2f_%2.2f_PixelMagnification_%2.3f_maxSF_%2.0fCPD.mat', stimPositionDegs(1), stimPositionDegs(2), rfMappingPixelMagnificationFactor, maxSFLimit);
    end

    coneMosaicSubspaceResponsesFileName = strrep(coneMosaicSubspaceResponsesFileName, '.mat', positionPostFix);
    mRGCMosaicSubspaceResponsesFileName = strrep(mRGCMosaicSubspaceResponsesFileName, '.mat', positionPostFix);
    optimallyMappedSubspaceRFmapsFileName = strrep(optimallyMappedSubspaceRFmapsFileName, '.mat', positionPostFix);

    if (reComputeInputConeMosaicSubspaceRFmappingResponses || reComputeMRGCMosaicSubspaceRFmappingResponses || reComputeRFs)

        % Compute the RF maps for ALL cells using stimuli at this 
        % position. Only some cells will be optimally mapped at each
        % position. These optimally derived RF maps are extracted by
        % the optimalyMappedRFsAtMosaicPosition() function
        computeRFmapsForAllCellsUsingStimuliAtTargetPosition(...
            theComputeReadyMRGCmosaic, opticsToEmploy, ...
            stimSizeDegs, stimPositionDegs, stimulusChromaticity, ...
            employConeFundamentalsDerivedFromInputConeMosaicAtStimPosition, ...
            maxSFLimit, rfMappingPixelMagnificationFactor, ...
            parPoolSize, ...
            coneMosaicSubspaceResponsesFileName, ...
            mRGCMosaicSubspaceResponsesFileName, ...
            optimallyMappedSubspaceRFmapsFileName, ...
            reComputeInputConeMosaicSubspaceRFmappingResponses, ...
            reComputeMRGCMosaicSubspaceRFmappingResponses, ...
            reComputeRFs, ...
            visualizedResponses);
    end

    if (visualizeOptimallyMappedRFmapLocations)
        visualizeAllOptimallyMappedRFmapLocations(optimallyMappedSubspaceRFmapsFileName, ...
            theComputeReadyMRGCmosaic, stimulusChromaticity, msequencePixelSizeDegs);
    end

end

function visualizeAllOptimallyMappedRFmapLocations(optimallyMappedSubspaceRFmapsFileName,...
    theMRGCMosaic, stimulusChromaticity, msequencePixelSizeDegs)
    % Save all the optimally mapped visual RF maps

    fprintf('Loading optimally mapped subspace RF maps from %s\n', optimallyMappedSubspaceRFmapsFileName);
    load(optimallyMappedSubspaceRFmapsFileName, 'optimallyMappedVisualRFmaps', 'indicesOfOptimallyMappedRGCs');


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

        [~,~,~, theCenterConeTypes, theCenterConeIndices] = ...
            theMRGCMosaic.centerConeTypeWeights(theRGCindex);

        theRGCRFcenterPositionDegs = theMRGCMosaic.rgcRFpositionsDegs(theRGCindex,:);

        [~, ~, theSurroundConeTypes, theSurroundConeIndices] = ...
            theMRGCMosaic.surroundConeTypeWeights(theRGCindex, theCenterConeIndices);

        idx = find(theSurroundConeTypes == cMosaic.LCONE_ID);
        surroundLconeIndices = theSurroundConeIndices(idx);

        idx = find(theSurroundConeTypes == cMosaic.MCONE_ID);
        surroundMconeIndices = theSurroundConeIndices(idx);

        d = optimallyMappedVisualRFmaps{iCell};
        d = generateMsequenceRFmap(d, msequencePixelSizeDegs);

        ax = subplot('Position', subplotPosVectors(1,1).v);
        renderSubspaceRFmap(ax, d, ...
            theMRGCMosaic.inputConeMosaic.coneRFpositionsDegs(surroundLconeIndices,:), ...
            theMRGCMosaic.inputConeMosaic.coneRFpositionsDegs(surroundMconeIndices,:), ...
            theRGCindex, theRGCRFcenterPositionDegs, stimulusChromaticity);

        xProfile = sum(d.theRFmap, 1);
        yProfile = sum(d.theRFmap, 2);
        maxP = max([max(abs(xProfile(:))) max(abs(yProfile(:)))]);

        ax = subplot('Position', subplotPosVectors(2,1).v);
        renderXProfile(ax, d.spatialSupportDegsX, xProfile, maxP, theRGCRFcenterPositionDegs);

        ax = subplot('Position', subplotPosVectors(1,2).v);
        renderYProfile(ax, d.spatialSupportDegsY, yProfile, maxP, theRGCRFcenterPositionDegs);

        ax = subplot('Position', subplotPosVectors(2,2).v);
        renderMsequenceRFmap(ax, d, theRGCRFcenterPositionDegs, msequencePixelSizeDegs);

        drawnow;
        pause
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

function renderXProfile(ax, spatialSupport, profile, maxP, theRGCRFcenterPositionDegs)
        cla(ax)

        shadedAreaPlotX(ax, spatialSupport, profile, 0, [0.5 0.5 0.5], [0.5 0.5 0.5], 0.5, 1.0);
        hold(ax, 'on');
        plot(ax, spatialSupport, profile*0, 'k-');

        plot(theRGCRFcenterPositionDegs(1)*[1 1], maxP*[-1 1], 'c-', 'LineWidth', 1.0);
        hold(ax, 'off');

        axis(ax, 'square')
        set(ax, 'XLim', [spatialSupport(1) spatialSupport(end)]);
        set(ax, 'YLim', maxP*[-1 1]);
        set(ax, 'XTick', -10:0.1:10);
        set(ax, 'XColor', [0 0 0], 'YColor', [1 0 0]);
        xlabel(ax, 'space, x (degs)');
        set(ax, 'FontSize', 16)
end

function renderYProfile(ax, spatialSupport, profile, maxP, theRGCRFcenterPositionDegs)
    cla(ax)
    yyaxis(ax, 'left');
    set(ax, 'YColor', 'none');
    yyaxis(ax, 'right');

    shadedAreaPlotY(ax, spatialSupport, profile, 0, [0.5 0.5 0.5], [0.5 0.5 0.5], 0.5, 1.0);
    hold(ax, 'on');
    plot(ax, profile*0, spatialSupport, 'k-');

    plot(maxP*[-1 1], theRGCRFcenterPositionDegs(2)*[1 1], 'c-', 'LineWidth', 1.0);
    hold(ax, 'off');

    axis(ax, 'square');
    set(ax, 'YLim', [spatialSupport(1) spatialSupport(end)]);
    set(ax, 'YTick', -10:0.1:10);
    set(ax, 'XDir', 'reverse', 'XLim', maxP*[-1 1]);
    set(ax, 'YColor', [0 0 0], 'XColor', [1 0 0]);
    ylabel(ax, 'space, y (degs)');
    set(ax, 'FontSize', 16);
end


function renderSubspaceRFmap(ax, d, surroundLconePositions, surroundMconePositions, theRGCindex, theRGCRFcenterPositionDegs, stimulusChromaticity)
    imagesc(ax, d.spatialSupportDegsX,  d.spatialSupportDegsY, d.theRFmap);
    hold(ax, 'on');    
    plot(ax, surroundLconePositions(:,1), surroundLconePositions(:,2), ...
         'r.', 'MarkerSize', 20, 'MarkerFaceColor', 'none', 'LineWidth', 1.0);
        
    plot(ax, surroundMconePositions(:,1), surroundMconePositions(:,2), ...
             'g.', 'MarkerSize', 20, 'MarkerFaceColor', 'none', 'LineWidth', 1.0);

    plot(theRGCRFcenterPositionDegs(1)*[1 1], [min(d.spatialSupportDegsY) max(d.spatialSupportDegsY)], 'c-', 'LineWidth', 1.0);
    plot([min(d.spatialSupportDegsX) max(d.spatialSupportDegsX)], theRGCRFcenterPositionDegs(2)*[1 1], 'c-', 'LineWidth', 1.0);
    hold(ax, 'off')
    axis(ax,'image'); axis(ax,'xy');
    set(ax, 'CLim', max(abs(d.theRFmap(:)))*[-1 1], ...
            'XLim', [d.spatialSupportDegsX(1) d.spatialSupportDegsX(end)], ...
            'YLim', [d.spatialSupportDegsY(1) d.spatialSupportDegsY(end)], ...
            'FontSize', 16 ...
    );
    set(ax, 'XTick', -10:0.1:10);
    set(ax, 'YTick', -10:0.1:10);
    colormap(ax,gray(1024));
    xlabel(ax, 'space, x (degs)');
    ylabel(ax, 'space, y (degs)');
    title(ax, sprintf('RGC %d - %s RF ([%2.2f ... %2.2f])', theRGCindex, stimulusChromaticity, min(d.theRFmap(:)), max(d.theRFmap(:))));
end


function renderMsequenceRFmap(ax, d, theRGCRFcenterPositionDegs, msequencePixelSizeDegs)
    imagesc(ax, d.theMsequenceRFmapSpatialSupportX, d.theMsequenceRFmapSpatialSupportY, d.theMsequenceRFmap);
    hold(ax, 'on');
    plot(theRGCRFcenterPositionDegs(1)*[1 1], [min(d.spatialSupportDegsY) max(d.spatialSupportDegsY)], 'c-', 'LineWidth', 1.0);
    plot([min(d.spatialSupportDegsX) max(d.spatialSupportDegsX)], theRGCRFcenterPositionDegs(2)*[1 1], 'c-', 'LineWidth', 1.0);
    hold(ax, 'off');
    axis(ax,'image'); axis(ax,'xy');
    set(ax, 'CLim', max(abs(d.theMsequenceRFmap(:)))*[-1 1], ...
            'XLim', [d.spatialSupportDegsX(1) d.spatialSupportDegsX(end)], ...
            'YLim', [d.spatialSupportDegsY(1) d.spatialSupportDegsY(end)], ...
            'FontSize', 16 ...
    );
    set(ax, 'XTick', -10:0.1:10);
    set(ax, 'YTick', -10:0.1:10);
    colormap(ax,gray(1024));
    xlabel(ax, 'space, x (degs)');
    ylabel(ax, 'space, y (degs)');
    title(ax, 'm-sequence mapped RF (pixel size: %2.3 arc min', msequencePixelSizeDegs*60));
end


function computeRFmapsForAllCellsUsingStimuliAtTargetPosition( ...
            theComputeReadyMRGCmosaic, opticsToEmploy, ...
            stimSizeDegs, stimPositionDegs, stimulusChromaticity,...
            employConeFundamentalsDerivedFromInputConeMosaicAtStimPosition, ...
            maxSFLimit, rfMappingPixelMagnificationFactor, ...
            parPoolSize,...
            coneMosaicSubspaceResponsesFileName, ...
            mRGCMosaicSubspaceResponsesFileName, ...
            optimallyMappedSubspaceRFmapsFileName, ...
            reComputeInputConeMosaicSubspaceRFmappingResponses, ...
            reComputeMRGCMosaicSubspaceRFmappingResponses, ...
            reComputeRFs, ...
            visualizedResponses)

    % Compute responses of the input cone mosaic to the subspace RF mapping stimuli
    if (reComputeInputConeMosaicSubspaceRFmappingResponses)

        fprintf('Cone mosaic subspace responses will be saved to %s \n', coneMosaicSubspaceResponsesFileName);
        MosaicPoolingOptimizer.generateInputConeMosaicSubspaceRFmappingLinearResponses(...
            theComputeReadyMRGCmosaic, opticsToEmploy, ...
            stimSizeDegs, stimPositionDegs, stimulusChromaticity, ...
            coneMosaicSubspaceResponsesFileName, ...
            'employConeFundamentalsDerivedFromInputConeMosaicAtStimPosition', employConeFundamentalsDerivedFromInputConeMosaicAtStimPosition, ...
            'maxSFLimit', maxSFLimit, ...
            'rfMappingPixelMagnificationFactor', rfMappingPixelMagnificationFactor, ...
            'visualizedResponses', visualizedResponses, ...
            'parPoolSize', parPoolSize);
    end

    fprintf('\nLoading cone mosaic subspace modulation responses and Hartley spatial modulation patterns ...');
    % Load the previously computed responses
    load(coneMosaicSubspaceResponsesFileName, ...
        'HartleySpatialModulationPatterns', 'spatialSupportDegs', 'stimParams', 'lIndices', 'mIndices', ...
        'theConeMosaicSubspaceLinearModulationResponses');
    HartleySpatialModulationPatterns = single(HartleySpatialModulationPatterns);
    fprintf('Done loading !\n');

    [HartleyStimNum, nCones] = size(theConeMosaicSubspaceLinearModulationResponses);

    if (reComputeMRGCMosaicSubspaceRFmappingResponses)
        fprintf('MRGC mosaic subspace RF maps and responses will be saved to %s \n', mRGCMosaicSubspaceResponsesFileName);
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

        parfor iStim = 1:HartleyStimNum
             fprintf('Computing mRGC mosaic response for Hartley pattern %d of %d (using %d parallel processes).\n', ...
                 iStim, HartleyStimNum, numWorkers);
             theConeMosaicModulationResponse = squeeze(theConeMosaicSubspaceLinearModulationResponses(iStim,:));
             theConeMosaicModulationResponse = reshape(theConeMosaicModulationResponse, [nTrials nTimePoints nCones]);

             % Compute mRGC mosaic responses based on the cone mosaic modulation responses
             [theMRGCMosaicResponse, ~, theMRGCresponseTemporalSupportSeconds] = ...
                    theComputeReadyMRGCmosaic.compute(theConeMosaicModulationResponse, theConeMosaicResponseTemporalSupportSeconds);
             theMRGCMosaicSubspaceRFmappingLinearResponses(iStim,:) = single(squeeze(theMRGCMosaicResponse(nTrials, nTimePoints ,:)));
        end
    
        fprintf('\nSaving computed mRGCRF mosaic SUBSPACE RF mapping linear responses to %s ...', mRGCMosaicSubspaceResponsesFileName);
        save(mRGCMosaicSubspaceResponsesFileName, ...
            'stimParams', ...
            'theMRGCMosaicSubspaceRFmappingLinearResponses', ...
            'spatialSupportDegs', 'lIndices', 'mIndices', ...
            '-v7.3');
    end

    if (reComputeRFs)
        % Load theMRGCMosaicSubspaceRFmappingLinearResponses
        load(mRGCMosaicSubspaceResponsesFileName, ...
                'theMRGCMosaicSubspaceRFmappingLinearResponses');
    
        indicesOfOptimallyMappedRGCs = indicesOfOptimallyMappedRGCsAtThisPosition(theComputeReadyMRGCmosaic, ...
            stimPositionDegs, stimSizeDegs);

        % Compute RF maps of all cells in the MRGC mosaic for stimuli delivered
        % in this position. Only some cells are optimally mapped in each
        % position grid, as the subspace stim size is usually small not covering the
        % entire RGC mosaic. These are selected by the
        % optimalyMappedRFsAtMosaicPosition() function later on
        theMRGCMosaicOptimallyMappedVisualRFmaps = computeRFs(...
            indicesOfOptimallyMappedRGCs, ...
            theMRGCMosaicSubspaceRFmappingLinearResponses, ...
            HartleySpatialModulationPatterns);
    
        fprintf('\nSaving computed visual RFs to %s ...', mRGCMosaicSubspaceResponsesFileName);
        save(mRGCMosaicSubspaceResponsesFileName, ...
            'theMRGCMosaicOptimallyMappedVisualRFmaps', 'indicesOfOptimallyMappedRGCs', '-append');
        fprintf('Done saving! \n');

        % Export visual RF maps for cells that are optimally mapped at this  position
        exportOptimalyMappedRFmaps(stimPositionDegs, mRGCMosaicSubspaceResponsesFileName, optimallyMappedSubspaceRFmapsFileName);
    end

end


function theRFmaps = computeRFs( ...
    indicesOfOptimallyMappedRGCs, ...
    theSubspaceRFmappingLinearResponses, ...
    HartleySpatialModulationPatterns)

    nStim = size(theSubspaceRFmappingLinearResponses,1);
    cellsNum = size(theSubspaceRFmappingLinearResponses,2);
    pixelsNum = size(HartleySpatialModulationPatterns,2);

    m = max(abs(theSubspaceRFmappingLinearResponses),[],1);
    cellsWithNonZeroResponse = find(m > 0);
    theRFmaps = cell(cellsNum, 1);

    parfor iCell = 1:numel(indicesOfOptimallyMappedRGCs)
        fprintf('Computing visual RF by accumulating Hartley patterns for the %d of %d optimally mapped RGC... \n', iCell, numel(indicesOfOptimallyMappedRGCs));
        theRGCindex = indicesOfOptimallyMappedRGCs(iCell);
        if (ismember(theRGCindex, cellsWithNonZeroResponse))
            theRFmap = zeros(pixelsNum, pixelsNum, 'single');
            allResponses = squeeze(theSubspaceRFmappingLinearResponses(:,theRGCindex));
            for iStim = 1:nStim
                theRFmap = theRFmap + ...
                        single(squeeze(HartleySpatialModulationPatterns(iStim,:,:)) * allResponses(iStim));
            end
            theRFmaps{iCell} = theRFmap;
        end
    end
end


function indicesOfOptimallyMappedRGCs = indicesOfOptimallyMappedRGCsAtThisPosition(theMRGCmosaic, ...
    stimPositionDegs, stimSizeDegs)

    % Find RGCs within the stimulus region
    if (numel(stimSizeDegs) == 2)
        widthDegs = stimSizeDegs(1);
        heightDegs = stimSizeDegs(2);
    else
        widthDegs = stimSizeDegs(1);
        heightDegs = stimSizeDegs(1);
    end

    theStimulusRegion = regionOfInterest(...
        'geometryStruct', struct(...
            'units', 'degs', ...
            'shape', 'rect', ...
            'center', stimPositionDegs, ...
            'width', widthDegs*0.75, ...
            'height', heightDegs*0.75 , ...
            'rotation', 0.0...
        ));

    indicesOfOptimallyMappedRGCs = theStimulusRegion.indicesOfPointsInside(theMRGCmosaic.rgcRFpositionsDegs);
    
    d = sqrt(sum(theMRGCmosaic.rgcRFpositionsDegs(indicesOfOptimallyMappedRGCs,:).^2,2));
    [~,idx] = sort(d, 'ascend');
    indicesOfOptimallyMappedRGCs = indicesOfOptimallyMappedRGCs(idx);

    % All the cells in the MRGC mosaic
    cellsNum = theMRGCmosaic.rgcsNum;
    
    fprintf('\nThere were %d/%d optimally mapped RF maps at (x,y) = %2.1f,%2.1f within the %2.2f x %2.2f mapping window\n', ...
        numel(indicesOfOptimallyMappedRGCs), cellsNum, ...
        stimPositionDegs(1), stimPositionDegs(2), widthDegs, heightDegs);

end

function exportOptimalyMappedRFmaps(stimPositionDegs, mRGCMosaicSubspaceResponsesFileName, optimallyMappedSubspaceRFmapsFileName)

    % Check to see if the subspace responses have been generated at this
    % grid position
    fileExists = isfile(mRGCMosaicSubspaceResponsesFileName);
    if (~fileExists)
        error('File %s does not exist.\n Not extracting optimally mapped visual RF maps at (x,y)=(%2.1f,%2.1f)\n', ...
            mRGCMosaicSubspaceResponsesFileName, ...
            stimPositionDegs(1), stimPositionDegs(2));
    end

    % Check to see if subspace RF maps have been computed at this grid position
    s = whos('-file', mRGCMosaicSubspaceResponsesFileName);
    computedVariableNames = cell(1, numel(s));
    for i = 1:numel(s)
        computedVariableNames{i} = s(i).name;
    end
    
    if (~ismember('theMRGCMosaicOptimallyMappedVisualRFmaps', computedVariableNames))
        error('Subspace responses are computed, however RF maps have not been computed yet.\nNot extracting optimally mapped visual RF maps at (x,y) = (%2.1f,%2.1f)\n', ...
            stimPositionDegs(1), stimPositionDegs(2));
    end
    
    fprintf('Loading RF map data from %s. Please wait ...\n', mRGCMosaicSubspaceResponsesFileName);
    % All good. Extract the computed subspace RF maps
    load(mRGCMosaicSubspaceResponsesFileName, ...
        'spatialSupportDegs', ...
        'indicesOfOptimallyMappedRGCs', ...
        'theMRGCMosaicOptimallyMappedVisualRFmaps');


    optimallyMappedRFsForThisGridPosition = numel(indicesOfOptimallyMappedRGCs);
    optimallyMappedVisualRFmaps = cell(1, optimallyMappedRFsForThisGridPosition);

    dx = spatialSupportDegs(2)-spatialSupportDegs(1);
    for iCell = 1:numel(theMRGCMosaicOptimallyMappedVisualRFmaps)
        optimallyMappedVisualRFmaps{iCell} = struct(...
            'theRFmap', theMRGCMosaicOptimallyMappedVisualRFmaps{iCell}, ...
            'spatialSupportDegsX', spatialSupportDegs+stimPositionDegs(1)-dx, ...
            'spatialSupportDegsY', spatialSupportDegs+stimPositionDegs(2)-dx ...
            );
    end


    % Save all the optimally mapped visual RF maps
    fprintf('Saving optimally mapped subspace RF maps to %s\n', optimallyMappedSubspaceRFmapsFileName);
    save(optimallyMappedSubspaceRFmapsFileName, ...
        'optimallyMappedVisualRFmaps', 'indicesOfOptimallyMappedRGCs', '-v7.3');

end