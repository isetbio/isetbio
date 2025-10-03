function computeVisualRFsUsingMSequenceMapping(...
        theComputeReadyMRGCmosaic, opticsToEmploy, ...
        stimSizeDegs, stimPositionDegs, rfPixelsAcross, ...
        ternaryInsteadOfBinaryMsequence, mSequenceBitLength, ...
        stimulusChromaticity, coneFundamentalsOptimizedForStimPosition, ...
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

    positionPostFix = sprintf('_atPosition_%2.2f_%2.2f.mat', stimPositionDegs(1), stimPositionDegs(2));

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
            ternaryInsteadOfBinaryMsequence, mSequenceBitLength, ...
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
            theComputeReadyMRGCmosaic, visualizedRGCindex, stimulusChromaticity, rfPixelsAcross);
    end

end

function computeRFmapsForAllCellsUsingStimuliAtTargetPosition( ...
            theComputeReadyMRGCmosaic, opticsToEmploy, ...
            stimSizeDegs, stimPositionDegs, rfPixelsAcross, ...
            stimulusChromaticity, coneFundamentalsOptimizedForStimPosition, ...
            ternaryInsteadOfBinaryMsequence, mSequenceBitLength, ...
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
            ternaryInsteadOfBinaryMsequence, mSequenceBitLength, ...
            stimulusChromaticity, coneFundamentalsOptimizedForStimPosition, ...
            coneMosaicResponsesFileName, ...
            'visualizedResponses', visualizedResponses, ...
            'parPoolSize', parPoolSize);
    end

    if (reComputeMRGCMosaicResponses)
        fprintf('\nLoading cone mosaic m-sequence modulation responses and spatial modulation patterns ...');
        % Load the previously computed responses
        load(coneMosaicResponsesFileName, ...
            'mSequenceIndicatorFunctions', 'spatialSupportDegs', 'stimParams', ...
            'theConeMosaicMSequenceLinearModulationResponses', ...
            'theConeMosaicMSequenceForwardModulationResponses', 'theConeMosaicMSequenceInverseModulationResponses');

        fprintf('Done loading !\n');

        [mSequenceStimsNum, nCones] = size(theConeMosaicMSequenceLinearModulationResponses);

        fprintf('MRGC mosaic m-sequence RF maps and responses will be saved to %s \n', mRGCMosaicResponsesFileName);
        fprintf('Computing visual m-sequence RF mapping responses for all RGCs in the mosaic ... \n');
        nTimePoints = 1;
        nTrials = 1;
        theConeMosaicResponseTemporalSupportSeconds = [0];
        
        theMRGCMosaicMSequenceRFmappingLinearResponses = zeros(mSequenceStimsNum, theComputeReadyMRGCmosaic.rgcsNum, 'single');
        theMRGCMosaicMSequenceRFmappingForwardResponses = zeros(mSequenceStimsNum, theComputeReadyMRGCmosaic.rgcsNum, 'single');
        theMRGCMosaicMSequenceRFmappingInverseResponses = zeros(mSequenceStimsNum, theComputeReadyMRGCmosaic.rgcsNum, 'single');

        % Use all processors
        [shutdownParPoolOnceCompleted, numWorkers] = MosaicPoolingOptimizer.resetParPool([]);
    
        % Noise-free responses
        theComputeReadyMRGCmosaic.noiseFlag = 'none';

        theMRGCresponseTemporalSupportSeconds = [];

        parfor iStim = 1:mSequenceStimsNum
             fprintf('Computing mRGC mosaic response for m-sequence frame %d of %d (using %d parallel processes).\n', ...
                 iStim, mSequenceStimsNum, numWorkers);

             % Compute mRGC mosaic linear responses 
             theConeMosaicModulationResponse = reshape(squeeze(theConeMosaicMSequenceLinearModulationResponses(iStim,:)), ...
                                                        [nTrials nTimePoints nCones]);
             [theMRGCMosaicResponse, ~, theMRGCresponseTemporalSupportSeconds] = ...
                    theComputeReadyMRGCmosaic.compute(theConeMosaicModulationResponse, theConeMosaicResponseTemporalSupportSeconds);
             theMRGCMosaicMSequenceRFmappingLinearResponses(iStim,:) = single(squeeze(theMRGCMosaicResponse(1,1,:)));

             % Repeat for forward polarity responses
             theConeMosaicModulationResponse = reshape(squeeze(theConeMosaicMSequenceForwardModulationResponses(iStim,:)), ...
                                                        [nTrials nTimePoints nCones]);
             theMRGCMosaicResponse = theComputeReadyMRGCmosaic.compute(theConeMosaicModulationResponse, theConeMosaicResponseTemporalSupportSeconds);
             theMRGCMosaicMSequenceRFmappingForwardResponses(iStim,:) = single(squeeze(theMRGCMosaicResponse(1,1,:)));

             % Repeat for inverse polarity responses
             theConeMosaicModulationResponse = reshape(squeeze(theConeMosaicMSequenceInverseModulationResponses(iStim,:)), ...
                                                        [nTrials nTimePoints nCones]);
             theMRGCMosaicResponse = theComputeReadyMRGCmosaic.compute(theConeMosaicModulationResponse, theConeMosaicResponseTemporalSupportSeconds);
             theMRGCMosaicMSequenceRFmappingInverseResponses(iStim,:) = single(squeeze(theMRGCMosaicResponse(1,1,:)));
        end
    
        fprintf('\nSaving computed mRGCRF mosaic SUBSPACE RF mapping linear responses to %s ...', mRGCMosaicResponsesFileName);
        save(mRGCMosaicResponsesFileName, ...
            'stimParams', ...
            'theMRGCMosaicMSequenceRFmappingLinearResponses', ...
            'theMRGCMosaicMSequenceRFmappingForwardResponses', ...
            'theMRGCMosaicMSequenceRFmappingInverseResponses', ...
            'theMRGCresponseTemporalSupportSeconds', ...
            'spatialSupportDegs', ...
            '-v7.3');
    end


    if (reComputeRFs)
        load(coneMosaicResponsesFileName, 'mSequenceIndicatorFunctions', 'spatialSupportDegs', 'stimParams');

        % Load theMRGCMosaicMSequenceRFmappingLinearResponses
        load(mRGCMosaicResponsesFileName, ...
            'theMRGCMosaicMSequenceRFmappingLinearResponses', ...
            'theMRGCMosaicMSequenceRFmappingForwardResponses', ...
            'theMRGCMosaicMSequenceRFmappingInverseResponses');
    
  
        % Determine indices of RGCs whose RF lie within the stimulus region
        indicesOfOptimallyMappedRGCs = MosaicPoolingOptimizer.indicesOfOptimallyMappedRGCsAtThisPosition(theComputeReadyMRGCmosaic, ...
            stimPositionDegs, stimSizeDegs);
        
        % Compute RF maps of all cells within the stimulus region
        [theMRGCMosaicOptimallyMappedVisualRFmaps, ...
         theMRGCMosaicOptimallyMappedVisualIncrementsRFmaps, ...
         theMRGCMosaicOptimallyMappedVisualDecrementsRFmaps] = computeRFs(...
            indicesOfOptimallyMappedRGCs, ...
            theMRGCMosaicMSequenceRFmappingLinearResponses, ...
            theMRGCMosaicMSequenceRFmappingForwardResponses, ...
            theMRGCMosaicMSequenceRFmappingInverseResponses, ...
            mSequenceIndicatorFunctions, ...
            stimParams.rfPixelRetinalPixelsWithin, ...
            spatialSupportDegs);
    
        fprintf('\nSaving computed visual RFs to %s ...', mRGCMosaicResponsesFileName);
        save(mRGCMosaicResponsesFileName, ...
            'theMRGCMosaicOptimallyMappedVisualRFmaps', ...
            'theMRGCMosaicOptimallyMappedVisualIncrementsRFmaps', ...
            'theMRGCMosaicOptimallyMappedVisualDecrementsRFmaps', ...
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
        'theMRGCMosaicOptimallyMappedVisualRFmaps', ...
        'theMRGCMosaicOptimallyMappedVisualIncrementsRFmaps', ...
        'theMRGCMosaicOptimallyMappedVisualDecrementsRFmaps');

    optimallyMappedVisualRFmaps = cell(1, numel(indicesOfOptimallyMappedRGCs));

    dx = spatialSupportDegs(2)-spatialSupportDegs(1);
    for iCell = 1:numel(theMRGCMosaicOptimallyMappedVisualRFmaps)
        optimallyMappedVisualRFmaps{iCell} = struct(...
            'theRFmap', theMRGCMosaicOptimallyMappedVisualRFmaps{iCell}, ...
            'theIncrementsRFmap', theMRGCMosaicOptimallyMappedVisualIncrementsRFmaps{iCell}, ...
            'theDecrementsRFmap', theMRGCMosaicOptimallyMappedVisualDecrementsRFmaps{iCell}, ...
            'spatialSupportDegsX', spatialSupportDegs+stimPositionDegs(1)-dx, ...
            'spatialSupportDegsY', spatialSupportDegs+stimPositionDegs(2)-dx);
    end

    % Save all the optimally mapped visual RF maps
    fprintf('Saving optimally mapped subspace RF maps to %s\n', optimallyMappedRFmapsFileName);
    save(optimallyMappedRFmapsFileName, ...
        'optimallyMappedVisualRFmaps', 'indicesOfOptimallyMappedRGCs', '-v7.3');

end

function [theRFmaps, theIncrementsRFmaps, theDecrementsRFmaps] = computeRFs(indicesOfOptimallyMappedRGCs, ...
    theLinearResponses, theForwardPolarityResponses, theInversePolarityResponses, ...
    mSequenceIndicatorFunctions, rfPixelRetinalPixelsWithin, spatialSupportDegs)

    nStim = size(theLinearResponses,1);
    cellsNum = size(theLinearResponses,2);
    pixelsNum = size(mSequenceIndicatorFunctions,2);

    m = max(abs(theLinearResponses),[],1);
    cellsWithNonZeroResponse = find(m > 0);
    theRFmaps = cell(cellsNum, 1);
    theIncrementsRFmaps = cell(cellsNum, 1);
    theDecrementsRFmaps = cell(cellsNum, 1);

    debugRFmap = ~true;

    parfor iCell = 1:numel(indicesOfOptimallyMappedRGCs)
        fprintf('Computing visual RF by accumulating m-sequence frames for the %d of %d optimally mapped RGC... \n', iCell, numel(indicesOfOptimallyMappedRGCs));
        theRGCindex = indicesOfOptimallyMappedRGCs(iCell);

        if (ismember(theRGCindex, cellsWithNonZeroResponse))
            theRFmap = zeros(pixelsNum, pixelsNum);
            theIncrementsRFmap = zeros(pixelsNum, pixelsNum);
            theDecrementsRFmap = zeros(pixelsNum, pixelsNum);

            allResponses = squeeze(theLinearResponses(:,theRGCindex));
            allForwardPolarityResponses = squeeze(theForwardPolarityResponses(:,theRGCindex));
            allInversePolarityResponses = squeeze(theInversePolarityResponses(:,theRGCindex));

            for iStim = 1:nStim
                theStimulusFrame = double(squeeze(mSequenceIndicatorFunctions(iStim,:,:)));
                theRFmap = theRFmap + theStimulusFrame * allResponses(iStim);
            end

            % Expand to full size
            theRFmaps{iCell} = rfMappingStimulusGenerator.expandFrame(theRFmap, rfPixelRetinalPixelsWithin);

            for iStim = 1:nStim
                theStimulusFrame = double(squeeze(mSequenceIndicatorFunctions(iStim,:,:)));

                theIncrementsStimulusFrame = theStimulusFrame*0;
                theIncrementsStimulusFrame(theStimulusFrame>0) = 1;

                theDecrementsStimulusFrame = theStimulusFrame*0;
                theDecrementsStimulusFrame(theStimulusFrame<0) = 1;

                theIncrementsRFmap = theIncrementsRFmap + ...
                    theIncrementsStimulusFrame * (allForwardPolarityResponses(iStim) - allInversePolarityResponses(iStim));

                theDecrementsRFmap = theDecrementsRFmap + ...
                    theDecrementsStimulusFrame * (allForwardPolarityResponses(iStim) - allInversePolarityResponses(iStim));

                if (debugRFmap)
                    m = max([max(abs(theIncrementsRFmap(:))) max(abs(theDecrementsRFmap(:)))]);

                    hFig = figure(44); clf;
                    set(hFig, 'Name', sprintf('%d/%d', iStim, nStim));
                    subplot(1,3,1)
                    imagesc(theIncrementsRFmap)
                    axis 'image'
                    set(gca, 'CLim', m*[-1 1]);
                    title(sprintf('max = %f', max(abs(theIncrementsRFmap(:)))));
                    colormap(gray)
        
                    subplot(1,3,2)
                    imagesc(theDecrementsRFmap)
                    axis 'image'
                    set(gca, 'CLim', m*[-1 1]);
                    title(sprintf('max = %f', max(abs(theDecrementsRFmap(:)))));
                    colormap(gray)
        
        
                    subplot(1,3,3)
                    compositeRFmap = theIncrementsRFmap + theDecrementsRFmap;
                    imagesc(compositeRFmap)
                    axis 'image'
                    set(gca, 'CLim', max(abs(compositeRFmap(:)))*[-1 1]);
                    title(sprintf('max = %f', max(abs(compositeRFmap(:)))));
                    colormap(gray)
        
                    drawnow
                end
            end
            
            % Expand to full size
            theIncrementsRFmaps{iCell} = rfMappingStimulusGenerator.expandFrame(theIncrementsRFmap, rfPixelRetinalPixelsWithin);
            theDecrementsRFmaps{iCell} = rfMappingStimulusGenerator.expandFrame(theDecrementsRFmap, rfPixelRetinalPixelsWithin);
        end

    end
end


function visualizeAllOptimallyMappedRFmapLocations(optimallyMappedRFmapsFileName,...
    theMRGCMosaic, visualizedRGCindex, stimulusChromaticity, rfPixelsAcross)

    % Load all the optimally mapped visual RF maps
    fprintf('Loading optimally mapped subspace RF maps from %s\n', optimallyMappedRFmapsFileName);
    load(optimallyMappedRFmapsFileName, 'optimallyMappedVisualRFmaps', 'indicesOfOptimallyMappedRGCs');

    hFig = figure(22); clf;
    set(hFig, 'Color', [1 1 1], 'Position', [10 10 1600 1110]);

    subplotPosVectors = NicePlot.getSubPlotPosVectors(...
       'rowsNum', 2, ...
       'colsNum', 3, ...
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

        [d.smoothedBipolarRF, smoothingKernel, rfPixelSizeSamples] = MosaicPoolingOptimizer.applyReidShapleySmoothingToRFmap(...
            d.spatialSupportDegsX, d.theRFmap, rfPixelsAcross);

        xLim = [d.spatialSupportDegsX(1) d.spatialSupportDegsX(end)];
        xLim = theRGCRetinalRFcenterPositionDegs(1) + 0.08*[-1 1];

        yLim = [d.spatialSupportDegsY(1) d.spatialSupportDegsY(end)];
        yLim = theRGCRetinalRFcenterPositionDegs(2) + 0.08*[-1 1];
        
        % The RF center position in the visual field
        [~,idx] = max(abs(d.smoothedBipolarRF(:)));
        [row,col] = ind2sub(size(d.theRFmap), idx);
        theRGCVisualRFcenterPositionDegs = [d.spatialSupportDegsX(col) d.spatialSupportDegsY(row)];

        crossHairsPosition = theRGCVisualRFcenterPositionDegs;
        %crossHairsPosition = theRGCRetinalRFcenterPositionDegs;

        ax = subplot('Position', subplotPosVectors(1,1).v);
        renderRFmap(ax, d, [], [], ...
            theMRGCMosaic.inputConeMosaic.coneRFpositionsDegs(theCenterConeIndices,:), ...
            theMRGCMosaic.inputConeMosaic.coneRFpositionsDegs(surroundLconeIndices,:), ...
            theMRGCMosaic.inputConeMosaic.coneRFpositionsDegs(surroundMconeIndices,:), ...
            theRGCindex, crossHairsPosition, ...
            stimulusChromaticity, xLim, yLim, ...
            'bipolar');

        ax = subplot('Position', subplotPosVectors(2,2).v);
        renderRFmap(ax, d, smoothingKernel, rfPixelSizeSamples, ...
            theMRGCMosaic.inputConeMosaic.coneRFpositionsDegs(theCenterConeIndices,:), ...
            theMRGCMosaic.inputConeMosaic.coneRFpositionsDegs(surroundLconeIndices,:), ...
            theMRGCMosaic.inputConeMosaic.coneRFpositionsDegs(surroundMconeIndices,:), ...
            theRGCindex, crossHairsPosition, ...
            stimulusChromaticity, xLim, yLim, ...
            'smoothedBipolarRF');


        ax = subplot('Position', subplotPosVectors(1,3).v);
        renderRFmap(ax, d, [], [],...
            theMRGCMosaic.inputConeMosaic.coneRFpositionsDegs(theCenterConeIndices,:), ...
            theMRGCMosaic.inputConeMosaic.coneRFpositionsDegs(surroundLconeIndices,:), ...
            theMRGCMosaic.inputConeMosaic.coneRFpositionsDegs(surroundMconeIndices,:), ...
            theRGCindex, crossHairsPosition, ...
            stimulusChromaticity, xLim, yLim, ...
            'increments');

        ax = subplot('Position', subplotPosVectors(2,3).v);
        renderRFmap(ax, d, [],[], ...
            theMRGCMosaic.inputConeMosaic.coneRFpositionsDegs(theCenterConeIndices,:), ...
            theMRGCMosaic.inputConeMosaic.coneRFpositionsDegs(surroundLconeIndices,:), ...
            theMRGCMosaic.inputConeMosaic.coneRFpositionsDegs(surroundMconeIndices,:), ...
            theRGCindex, crossHairsPosition, ...
            stimulusChromaticity, xLim, yLim, ...
            'decrements');


        xProfile = sum(d.theRFmap, 1);
        xIncrementsProfile = sum(d.theIncrementsRFmap, 1);
        xDecrementsProfile = sum(d.theDecrementsRFmap, 1);

        yProfile = sum(d.theRFmap, 2);
        yIncrementsProfile = sum(d.theIncrementsRFmap, 2);
        yDecrementsProfile = sum(d.theDecrementsRFmap, 2);

        maxP = max([...
            max(abs(xProfile(:))) ...
            max(abs(xIncrementsProfile(:))) ...
            max(abs(xDecrementsProfile(:))) ...
            max(abs(yProfile(:))) ...
            max(abs(yIncrementsProfile(:))) ...
            max(abs(yDecrementsProfile(:))) ...
            ]);

        ax = subplot('Position', subplotPosVectors(2,1).v);
        renderXProfile(ax, d.spatialSupportDegsX, xProfile, xIncrementsProfile, xDecrementsProfile, maxP, crossHairsPosition, xLim);

        ax = subplot('Position', subplotPosVectors(1,2).v);
        renderYProfile(ax, d.spatialSupportDegsY, yProfile, yIncrementsProfile, yDecrementsProfile, maxP, crossHairsPosition, yLim);

        drawnow;
        disp('Hit enter to continue : ')
        pause;
    end
end


function renderXProfile(ax, spatialSupport, profile, incrementsProfile, decrementsProfile, maxP, crossHairsPosition, xLim)

    cla(ax)

    shadedAreaPlotX(ax, spatialSupport, profile, 0, [0.5 0.5 0.5], [0.5 0.5 0.5], 0.5, 1.0);
    hold(ax, 'on');
    plot(ax, spatialSupport, profile*0, 'k-', 'LineWidth', 1.0);
    plot(ax, spatialSupport, incrementsProfile, 'r-', 'LineWidth', 1.5);
    plot(ax, spatialSupport, decrementsProfile, 'b-', 'LineWidth', 1.5);

    plot(crossHairsPosition(1)*[1 1], maxP*[-1 1], 'b-', 'LineWidth', 1.0);
    hold(ax, 'off');

    axis(ax, 'square')
    set(ax, 'XLim', xLim);
    set(ax, 'YLim', maxP*[-1 1]);
    set(ax, 'XTick', -10:0.05:10);
    set(ax, 'XColor', [0 0 0], 'YColor', [1 0 0]);
    xlabel(ax, 'space, x (degs)');
    set(ax, 'FontSize', 16)
end

function renderYProfile(ax, spatialSupport, profile, incrementsProfile, decrementsProfile, maxP, crossHairsPosition, yLim)
    cla(ax)
    yyaxis(ax, 'left');
    set(ax, 'YColor', 'none');
    yyaxis(ax, 'right');

    shadedAreaPlotY(ax, spatialSupport, profile, 0, [0.5 0.5 0.5], [0.5 0.5 0.5], 0.5, 1.0);
    hold(ax, 'on');
    plot(ax, profile*0, spatialSupport, 'k-', 'LineWidth', 1.0);
    plot(ax, incrementsProfile, spatialSupport, 'r-', 'LineWidth', 1.5);
    plot(ax, decrementsProfile, spatialSupport, 'b-', 'LineWidth', 1.5);

    plot(maxP*[-1 1], crossHairsPosition(2)*[1 1], 'b-', 'LineWidth', 1.0);
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


function renderRFmap(ax, d, smoothingKernel, rfPixelSizeSamples, centerConePositions, surroundLconePositions, surroundMconePositions, theRGCindex, ...
    crossHairsPosition, stimulusChromaticity, xLim, yLim,  rfPolarity)

    switch (rfPolarity)
        case 'smoothedBipolarRF'
            theRFmap = d.smoothedBipolarRF;
            rfMapRange = max(abs(d.theRFmap(:)))*[-1 1];

        case 'bipolar'
            theRFmap = d.theRFmap;
            rfMapRange = max(abs(d.theRFmap(:)))*[-1 1];

        case 'increments'
            theRFmap = d.theIncrementsRFmap;
            rfMapRange = max([ max(abs(d.theIncrementsRFmap(:))) max(abs(d.theDecrementsRFmap(:)))])*[-1 1];

        case 'decrements'
            theRFmap = d.theDecrementsRFmap;
            rfMapRange = max([ max(abs(d.theIncrementsRFmap(:))) max(abs(d.theDecrementsRFmap(:)))])*[-1 1];

        otherwise
            error('Unknown RF polarity: ''%s''.', rfPolarity);
    end

    imagesc(ax, d.spatialSupportDegsX,  d.spatialSupportDegsY, theRFmap);
    hold(ax, 'on');    

    dx = d.spatialSupportDegsX(2)- d.spatialSupportDegsX(1);

    if (~isempty(smoothingKernel))
        imagesc(ax, xLim(1) + dx*(1:size(smoothingKernel,2)), ...
                    yLim(1)+dx*(1:size(smoothingKernel,1)), ...
                    rfMapRange(2)*smoothingKernel/max(smoothingKernel(:)));
    end

    if (~isempty(rfPixelSizeSamples))
        
        
        xx(1) = mean(d.spatialSupportDegsX);
        k = 0;
        while (xx(1)+k*(rfPixelSizeSamples)*dx<d.spatialSupportDegsX(end))
            xx(numel(xx)+1) = xx(1) + k*(rfPixelSizeSamples)*dx;
            xx(numel(xx)+1) = xx(1) - k*(rfPixelSizeSamples)*dx;
            k = k + 1;
        end
        
        yy(1) = mean(d.spatialSupportDegsY);
        k = 0;
        while (yy(1)+k*(rfPixelSizeSamples)*dx<d.spatialSupportDegsY(end))
            yy(numel(yy)+1) = yy(1) + k*(rfPixelSizeSamples)*dx;
            yy(numel(yy)+1) = yy(1) - k*(rfPixelSizeSamples)*dx;
            k = k + 1;
        end

        for i = 1:numel(xx)
            plot(ax, xx(i)*[1 1], [d.spatialSupportDegsY(1) d.spatialSupportDegsY(end)], 'k-');
        end
        for i = 1:numel(yy)
            plot(ax,  [d.spatialSupportDegsX(1) d.spatialSupportDegsX(end)], yy(i)*[1 1], 'k-');
        end

    else
        scatter(ax, surroundLconePositions(:,1), surroundLconePositions(:,2), ...
             16*16, 'MarkerFaceColor', [1 0.5 0.5], 'MarkerEdgeColor', [1 0 0], 'MarkerFaceAlpha', 0.0, 'LineWidth', 1.0);
            
        scatter(ax, surroundMconePositions(:,1), surroundMconePositions(:,2), ...
             16*16, 'MarkerFaceColor', [0.5 0.8 0.5], 'MarkerEdgeColor', [0 0.8 0], 'MarkerFaceAlpha', 0.0, 'LineWidth', 1.0);
    end

    plot(ax, mean(centerConePositions(:,1),1), mean(centerConePositions(:,2),1), 'kx', 'MarkerSize', 16, 'LineWidth', 2.0);

    plot(crossHairsPosition(1)*[1 1], [-100 100], 'b-', 'LineWidth', 1.0);
    plot([-100 100], crossHairsPosition(2)*[1 1], 'b-', 'LineWidth', 1.0);

    hold(ax, 'off')
    axis(ax,'image'); axis(ax,'xy');
    set(ax, 'CLim', rfMapRange, ...
            'XLim', xLim, ...
            'YLim', yLim, ...
            'FontSize', 16 ...
    );
    set(ax, 'XTick', -10:0.05:10);
    set(ax, 'YTick', -10:0.05:10);
    
    cLUT = MosaicPoolingOptimizer.generateReidShapleyRFmapLUT();

    colormap(ax, cLUT);
    midPoint = (size(cLUT,1)-1)/2+1;
    colorbar
    set(ax, 'Color', cLUT(midPoint,:));
    xlabel(ax, 'space, x (degs)');
    ylabel(ax, 'space, y (degs)');
    title(ax, sprintf('RGC %d - %s %s RF ([%2.2f ... %2.2f])', theRGCindex, stimulusChromaticity, rfPolarity, min(theRFmap(:)), max(theRFmap(:))));
end


