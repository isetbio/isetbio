function computeVisualRFsUsingMSequenceMapping(...
        theComputeReadyMRGCmosaic, opticsToEmploy, ...
        stimSizeDegs, stimPositionDegs, rfPixelsAcross, ...
        stimulusChromaticity, coneFundamentalsOptimizedForStimPosition, ...
        ternaryInsteadOfBinaryMsequence, ...
        coneMosaicResponsesFileName, ...
        mRGCMosaicResponsesFileName, ...
        optimallyMappedRFmapsFileName, ...
        reComputeInputConeMosaicResponses, ...
        reComputeMRGCMosaicResponses, ...
        reComputeRFs, ...
        visualizeOptimallyMappedRFmapLocations, varargin)

    p = inputParser;
    p.addParameter('parPoolSize', [], @(x)(isempty(x)||(isscalar(x))));
    p.addParameter('visualizedResponses', false, @islogical);
    p.addParameter('visualizedRGCindex', [], @(x)(isempty(x)||(isscalar(x))));
    p.parse(varargin{:});

    parPoolSize = p.Results.parPoolSize;
    visualizedRGCindex = p.Results.visualizedRGCindex;
    visualizedResponses = p.Results.visualizedResponses;

    if (isempty(stimPositionDegs))
        stimPositionDegs = theComputeReadyMRGCmosaic.eccentricityDegs;
    end

    positionPostFix = sprintf('_atPosition_%2.2f_%2.2f_PixelsAcross_%d.mat', stimPositionDegs(1), stimPositionDegs(2), rfPixelsAcross);

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
            stimSizeDegs, stimPositionDegs, rfPixelsAcross, ...
            stimulusChromaticity, coneFundamentalsOptimizedForStimPosition, ...
            ternaryInsteadOfBinaryMsequence, ...
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
            theComputeReadyMRGCmosaic, visualizedRGCindex, stimulusChromaticity);
    end

end

function computeRFmapsForAllCellsUsingStimuliAtTargetPosition( ...
            theComputeReadyMRGCmosaic, opticsToEmploy, ...
            stimSizeDegs, stimPositionDegs, rfPixelsAcross, ...
            stimulusChromaticity, coneFundamentalsOptimizedForStimPosition, ...
            ternaryInsteadOfBinaryMsequence, ...
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
        fprintf('Cone mosaic m-sequence responses will be saved to %s \n', coneMosaicResponsesFileName);
        MosaicPoolingOptimizer.generateInputConeMosaicMSequenceRFmappingLinearResponses(...
            theComputeReadyMRGCmosaic, opticsToEmploy, ...
            stimSizeDegs, stimPositionDegs, rfPixelsAcross, ...
            stimulusChromaticity, coneFundamentalsOptimizedForStimPosition, ...
            ternaryInsteadOfBinaryMsequence, ...
            coneMosaicResponsesFileName, ...
            'visualizedResponses', visualizedResponses, ...
            'parPoolSize', parPoolSize);
    end

    if (reComputeMRGCMosaicResponses)
        fprintf('\nLoading cone mosaic m-sequence modulation responses and spatial modulation patterns ...');
        % Load the previously computed responses
        load(coneMosaicResponsesFileName, ...
            'mSequenceSpatialModulationPatterns', 'spatialSupportDegs', 'stimParams', ...
            'theConeMosaicMSequenceLinearModulationResponses');
        mSequenceSpatialModulationPatterns = single(mSequenceSpatialModulationPatterns);
        fprintf('Done loading !\n');

        [mSequenceStimsNum, nCones] = size(theConeMosaicMSequenceLinearModulationResponses);

        fprintf('MRGC mosaic m-sequence RF maps and responses will be saved to %s \n', mRGCMosaicResponsesFileName);
        fprintf('Computing visual subspace RF mapping responses for all RGCs in the mosaic ... \n');
        nTimePoints = 1;
        nTrials = 1;
        theConeMosaicResponseTemporalSupportSeconds = [0];
        
        theMRGCMosaicMSequenceRFmappingLinearResponses = zeros(...
                     mSequenceStimsNum, theComputeReadyMRGCmosaic.rgcsNum, 'single');

        % Use all processors
        [shutdownParPoolOnceCompleted, numWorkers] = MosaicPoolingOptimizer.resetParPool([]);
    
        % Noise-free responses
        theComputeReadyMRGCmosaic.noiseFlag = 'none';

        parfor iStim = 1:mSequenceStimsNum
             fprintf('Computing mRGC mosaic response for m-sequence frame %d of %d (using %d parallel processes).\n', ...
                 iStim, mSequenceStimsNum, numWorkers);
             theConeMosaicModulationResponse = squeeze(theConeMosaicMSequenceLinearModulationResponses(iStim,:));
             theConeMosaicModulationResponse = reshape(theConeMosaicModulationResponse, [nTrials nTimePoints nCones]);

             % Compute mRGC mosaic responses based on the cone mosaic modulation responses
             [theMRGCMosaicResponse, ~, theMRGCresponseTemporalSupportSeconds] = ...
                    theComputeReadyMRGCmosaic.compute(theConeMosaicModulationResponse, theConeMosaicResponseTemporalSupportSeconds);
             theMRGCMosaicMSequenceRFmappingLinearResponses(iStim,:) = single(squeeze(theMRGCMosaicResponse(nTrials, nTimePoints ,:)));
        end
    
        fprintf('\nSaving computed mRGCRF mosaic SUBSPACE RF mapping linear responses to %s ...', mRGCMosaicResponsesFileName);
        save(mRGCMosaicResponsesFileName, ...
            'stimParams', ...
            'theMRGCMosaicMSequenceRFmappingLinearResponses', ...
            'spatialSupportDegs', ...
            '-v7.3');
    end


    if (reComputeRFs)
        load(coneMosaicResponsesFileName, 'mSequenceSpatialModulationPatterns', 'spatialSupportDegs', 'stimParams');
        mSequenceSpatialModulationPatterns = single(mSequenceSpatialModulationPatterns);

        % Load theMRGCMosaicMSequenceRFmappingLinearResponses
        load(mRGCMosaicResponsesFileName, 'theMRGCMosaicMSequenceRFmappingLinearResponses');
    
  
        % Determine indices of RGCs whose RF lie within the stimulus region
        indicesOfOptimallyMappedRGCs = MosaicPoolingOptimizer.indicesOfOptimallyMappedRGCsAtThisPosition(theComputeReadyMRGCmosaic, ...
            stimPositionDegs, stimSizeDegs);
        
        % Compute RF maps of all cells within the stimulus region
        theMRGCMosaicOptimallyMappedVisualRFmaps = computeRFs(...
            indicesOfOptimallyMappedRGCs, ...
            theMRGCMosaicMSequenceRFmappingLinearResponses, ...
            mSequenceSpatialModulationPatterns, ...
            spatialSupportDegs);
    
        fprintf('\nSaving computed visual RFs to %s ...', mRGCMosaicResponsesFileName);
        save(mRGCMosaicResponsesFileName, ...
            'theMRGCMosaicOptimallyMappedVisualRFmaps', ...
            'indicesOfOptimallyMappedRGCs', '-append');
        fprintf('Done saving! \n');

        % Export visual RF maps for cells that are optimally mapped at this  position
        exportOptimalyMappedRFmaps(stimPositionDegs, mRGCMosaicResponsesFileName, optimallyMappedRFmapsFileName);
    end

end

function exportOptimalyMappedRFmaps(stimPositionDegs, mRGCMosaicResponsesFileName, optimallyMappedRFmapsFileName)
    
    fprintf('Loading RF map data from %s. Please wait ...\n', mRGCMosaicResponsesFileName);
    % All good. Extract the computed subspace RF maps
    load(mRGCMosaicResponsesFileName, ...
        'spatialSupportDegs', ...
        'indicesOfOptimallyMappedRGCs', ...
        'theMRGCMosaicOptimallyMappedVisualRFmaps');

    optimallyMappedVisualRFmaps = cell(1, numel(indicesOfOptimallyMappedRGCs));

    dx = spatialSupportDegs(2)-spatialSupportDegs(1);
    for iCell = 1:numel(theMRGCMosaicOptimallyMappedVisualRFmaps)
        optimallyMappedVisualRFmaps{iCell} = struct(...
            'theRFmap', theMRGCMosaicOptimallyMappedVisualRFmaps{iCell}, ...
            'spatialSupportDegsX', spatialSupportDegs+stimPositionDegs(1)-dx, ...
            'spatialSupportDegsY', spatialSupportDegs+stimPositionDegs(2)-dx);
    end

    % Save all the optimally mapped visual RF maps
    fprintf('Saving optimally mapped subspace RF maps to %s\n', optimallyMappedRFmapsFileName);
    save(optimallyMappedRFmapsFileName, ...
        'optimallyMappedVisualRFmaps', 'indicesOfOptimallyMappedRGCs', '-v7.3');

end

function theRFmaps = computeRFs(indicesOfOptimallyMappedRGCs, ...
    theRFmappingLinearResponses, spatialModulationPatterns, spatialSupportDegs)

    nStim = size(theRFmappingLinearResponses,1);
    cellsNum = size(theRFmappingLinearResponses,2);
    pixelsNum = size(spatialModulationPatterns,2);

    m = max(abs(theRFmappingLinearResponses),[],1);
    cellsWithNonZeroResponse = find(m > 0);
    theRFmaps = cell(cellsNum, 1);

    parfor iCell = 1:numel(indicesOfOptimallyMappedRGCs)
        fprintf('Computing visual RF by accumulating m-sequence frames for the %d of %d optimally mapped RGC... \n', iCell, numel(indicesOfOptimallyMappedRGCs));
        theRGCindex = indicesOfOptimallyMappedRGCs(iCell);

        if (ismember(theRGCindex, cellsWithNonZeroResponse))
            theRFmap = zeros(pixelsNum, pixelsNum, 'single');
            allResponses = squeeze(theRFmappingLinearResponses(:,theRGCindex));

            for iStim = 1:nStim
                theRFmap = theRFmap + single(squeeze(spatialModulationPatterns(iStim,:,:)) * allResponses(iStim));
            end
            theRFmaps{iCell} = theRFmap;
        end

    end
end


function visualizeAllOptimallyMappedRFmapLocations(optimallyMappedRFmapsFileName,...
    theMRGCMosaic, visualizedRGCindex, stimulusChromaticity)

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

        xLim = [d.spatialSupportDegsX(1) d.spatialSupportDegsX(end)];
        xLim = theRGCRetinalRFcenterPositionDegs(1) + 0.08*[-1 1];

        yLim = [d.spatialSupportDegsY(1) d.spatialSupportDegsY(end)];
        yLim = theRGCRetinalRFcenterPositionDegs(2) + 0.08*[-1 1];

        % The RF center position in the visual field
        [~,idx] = max(abs(d.theRFmap(:)));
        [row,col] = ind2sub(size(d.theRFmap), idx);
        theRGCVisualRFcenterPositionDegs = [d.spatialSupportDegsX(col) d.spatialSupportDegsY(row)];

        ax = subplot('Position', subplotPosVectors(1,1).v);
        renderRFmap(ax, d, ...
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
        renderYProfile(ax, d.spatialSupportDegsY, yProfile, maxP, theRGCVisualRFcenterPositionDegs, yLim);

        drawnow;
    end
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


function renderRFmap(ax, d, centerConePositions, surroundLconePositions, surroundMconePositions, theRGCindex, ...
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


