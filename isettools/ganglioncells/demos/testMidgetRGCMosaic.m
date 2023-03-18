function testMidgetRGCMosaic
    
    arbitraryNodesToCompute =  selectNodesToRecompute();
    mosaicEcc = 2.5;
    mosaicEcc = 7.0;
    %mosaicEcc = -10.0;

    % Get mosaic ecc and size
    mosaicParams = getMosaicParams(mosaicEcc);

    % Ask user what operation to perform
    operationSetToPerformContains = promptUserForOperationsToPerform(mosaicParams);
    
   
    % Generate mosaic filename
    [mosaicFileName, resourcesDirectory] = ...
        MosaicPoolingOptimizer.resourceFileNameAndPath('mosaic', ...
            'mosaicParams', mosaicParams);

    if (operationSetToPerformContains.generateRGCMosaic)
        % Generate mosaic, its input coneMosaic and connect cones to the RF centers
        theMidgetRGCMosaic = mRGCMosaic(...
            'eccentricityDegs', mosaicParams.eccDegs, ...
            'sizeDegs', mosaicParams.sizeDegs);
    
        % Save the generated center-only connected mRGCmosaic
        save(fullfile(resourcesDirectory,mosaicFileName), 'theMidgetRGCMosaic');

        % Visualize input cone mosaic
        theMidgetRGCMosaic.inputConeMosaic.visualize()
    
        % Visualize cone pooling by the RF centers
        theMidgetRGCMosaic.visualize();

    else
        % Load a previously-generated mRGCMosaic
        load(fullfile(resourcesDirectory,mosaicFileName), 'theMidgetRGCMosaic');
    end



    if (operationSetToPerformContains.visualizeRGCMosaic)
        % Load the generated center-only connected mRGCmosaic
        load(fullfile(resourcesDirectory,mosaicFileName), 'theMidgetRGCMosaic');

        identifyPooledCones = true;
        identifyInputCones = true;
        plotRFoutlines = true;
    
        hFig = figure(1); clf;
        set(hFig, 'Position', [10 10 1800 1000], 'Color', [1 1 1]);
        theMidgetRGCMosaic.visualize(...
            'figureHandle', hFig, ...
            'identifiedConeApertureThetaSamples', 32, ...
            'identifiedConeAperture', 'lightCollectingAreaCharacteristicDiameter', ...
            'identifyPooledCones', identifyPooledCones, ...
            'identifyInputCones',identifyInputCones, ...
            'plotRFoutlines', plotRFoutlines, ...
            'domainVisualizationLimits', [], ...
            'domainVisualizationTicks', [], ...
            'backgroundColor', [0 0 0]);

        % Report stats
        centerConesNumCases = theMidgetRGCMosaic.centerConePoolingStats();

        % Ask user whether to eliminate certain RGCs
        eliminateRGCs = input('Eliminate RGCs with a certain no of center cones? [y = YES] : ', 's');
        if (strcmpi(eliminateRGCs, 'y'))
            % Ask user which RGCs to eliminate
            targetCenterConesNum = input(sprintf('What number of center cones [%d .. %d]? :', min(centerConesNumCases), max(centerConesNumCases)));
            
            % Eliminate the target RGCs
            theMidgetRGCMosaic.eliminateRGCsWithThisManyCenterConesNum(targetCenterConesNum);
            
            % Ask user whether to save the updated RGC mosaic
            exportUpdatedRGCMosaic = input('Export the updated RGC mosaic [y = YES] : ', 's');
            if (strcmpi(exportUpdatedRGCMosaic, 'y'))
                save(fullfile(resourcesDirectory,mosaicFileName), 'theMidgetRGCMosaic', '-v7.3');
                fprintf('The updated RGC mosaic was saved to %s.\n', ...
                    fullfile(resourcesDirectory,mosaicFileName));
            else
                fprintf('The updated RGC mosaic was not saved to the disk.\n');
            end
        else
            fprintf('No change in the RGC mosaic\n');
        end

    end


    % Default optics params
    opticsParams = theMidgetRGCMosaic.defaultOpticsParams;

    if (operationSetToPerformContains.visualizePSFsAtEccentricities)
        xPosDegs = 0.5*mosaicParams.sizeDegs(1)*[-1 0 1];
        yPosDegs = 0.5*mosaicParams.sizeDegs(2)*[-1 0 1];

        [X,Y] = meshgrid(xPosDegs, yPosDegs);
        eccDegs = [X(:) Y(:)];
        eccDegs = bsxfun(@plus, eccDegs, mosaicParams.eccDegs);

        theMidgetRGCMosaic.visualizeOpticsAtEccentricities(eccDegs, opticsParams);
    end


    % RetinalRFmodel params to employ
    % Change something if we want, like the model name, e.g. choose cell index 3,
    % 'arbitraryCenterConeWeights_doubleExpH1cellIndex3SurroundWeights', ... 
    H1cellIndex = 4;
    retinalRFmodelParams = MosaicPoolingOptimizer.defaultRetinalRFmodelParams;
    retinalRFmodelParams.conePoolingModel = sprintf('arbitraryCenterConeWeights_doubleExpH1cellIndex%dSurroundWeights', H1cellIndex);

    % Generate filename for the computed coneMosaicSTF responses
    [coneMosaicSTFresponsesFileName, resourcesDirectory] = ...
        MosaicPoolingOptimizer.resourceFileNameAndPath('coneMosaicSTFresponses', ...
            'mosaicParams', mosaicParams, ...
            'opticsParams', opticsParams);
   
    % Generate filename for the computed mRGCMosaicSTF responses
    [mRGCMosaicSTFresponsesFileName, resourcesDirectory] = ...
        MosaicPoolingOptimizer.resourceFileNameAndPath('mRGCMosaicSTFresponses', ...
            'mosaicParams', mosaicParams, ...
            'opticsParams', opticsParams);

    % Generate filename for the computed coneMosaic subspace responses
    [coneMosaicSubspaceResponsesFileName, resourcesDirectory] = ...
        MosaicPoolingOptimizer.resourceFileNameAndPath('coneMosaicSubspaceResponses', ...
            'mosaicParams', mosaicParams, ...
            'opticsParams', opticsParams);

    % Generate filename for the computed mRGCMosaic subspace responses
    [mRGCMosaicSubspaceRresponsesFileName, resourcesDirectory] = ...
        MosaicPoolingOptimizer.resourceFileNameAndPath('mRGCMosaicSubspaceResponses', ...
            'mosaicParams', mosaicParams, ...
            'opticsParams', opticsParams);


    % Generate the filenames of the optimizedRGCpooling objects
    [optimizedRGCpoolingObjectsFileName, resourcesDirectory] = ...
            MosaicPoolingOptimizer.resourceFileNameAndPath('optimizedRGCpoolingObjects', ...
                'mosaicParams', mosaicParams, ...
                'opticsParams', opticsParams, ...
                'retinalRFmodelParams', retinalRFmodelParams);

    % Generate the filename of the compute-ready mRGCMosaic to generate
    [computeReadyMosaicFileName, resourcesDirectory] = ...
            MosaicPoolingOptimizer.resourceFileNameAndPath('computeReadyMosaic', ...
                'mosaicParams', mosaicParams, ...
                'opticsParams', opticsParams, ...
                'retinalRFmodelParams', retinalRFmodelParams);


    % Stage 2: Generate optics and compute input cone mosaic STF responses
    if (operationSetToPerformContains.computeConeMosaicSTFresponses)

        opticsChoice = 'f';
        validOpticsChoices = {'n', 'c'};
        while (~ismember(opticsChoice, validOpticsChoices))
            fprintf('\nOptics options: \n');
            fprintf('[n] - Native optics (at the mosaic''s center) \n');
            fprintf('[c] - Custom optics (at an arbitrary position) \n');
            opticsChoice = input('Optics to employ: ', 's');
        end

        switch (opticsChoice)
            case 'n'
                % Generate the native optics appropriate for the midgetRGCMosaic at hand. 
                % These will be used to compute the
                % cone mosaic STF responses and optimize the RGC RFs
                theMidgetRGCMosaic.generateNativeOptics(opticsParams);
                opticsToEmploy = 'native';
            case 'c'
                xyPos = input('Enter (x,y) position (in degrees) for the optics (e.g., [2 0]) : ');
                opticsPositionPostfix = sprintf('AtCustomXYposition_%2.2f_%2.2f.mat', xyPos(1), xyPos(2));
                coneMosaicSTFresponsesFileName = strrep(coneMosaicSTFresponsesFileName, '.mat', opticsPositionPostfix);
                
                % Generate optics at the user-specified position
                opticsParams.positionDegs = xyPos;
                theMidgetRGCMosaic.generateNativeOptics(opticsParams);
                opticsToEmploy = 'custom';

            otherwise
                fprintf('Valid optics options: ''n'' (native, at mosaic''s center), or ''c'', (custom position)\n');
                error('Unknown optics choice: ''%s''.', opticsChoice)
        end

        compuseConeMosaicSTFResponsesSeparatelyAtEachGridLocation = false;
        if (compuseConeMosaicSTFResponsesSeparatelyAtEachGridLocation)
            % Instantiate the mosaic pooling optimizer
            theMosaicPoolingOptimizer = MosaicPoolingOptimizer(...
                theMidgetRGCMosaic, ...
                'generateSamplingGrids', true);

            % Compute the input cone mosaic visual STFs, at a single
            % gridNode (using small stimulus patches centered on that node)
            for gridNodeIndex = 1:theMosaicPoolingOptimizer.gridNodesNum
                % Responses filename
                coneMosaicSTFresponsesFileNameForThisNode = strrep(coneMosaicSTFresponsesFileName, '.mat', sprintf('_AtGridNode_%d.mat', gridNodeIndex));
    
                stimSizeDegs = [1 1];
                theMosaicPoolingOptimizer.generateInputConeMosaicSTFresponses(...
                    gridNodeIndex, stimSizeDegs,  ...
                    fullfile(resourcesDirectory, coneMosaicSTFresponsesFileNameForThisNode), ...
                    'useParfor', true, ...
                    'visualizedResponses', ~true, ...
                    'opticsToEmploy', opticsToEmploy);
            end
        else
            % Instantiate the mosaic pooling optimizer
            theMosaicPoolingOptimizer = MosaicPoolingOptimizer(...
                theMidgetRGCMosaic, ...
                'generateSamplingGrids', false);

            % Compute the input cone mosaic visual STFs across all nodes, 
            % using full field stimuli
            theMosaicPoolingOptimizer.generateInputConeMosaicSTFresponses(...
                    [], [], ...
                    fullfile(resourcesDirectory, coneMosaicSTFresponsesFileName), ...
                    'useParfor', ~true, ...
                    'visualizedResponses', ~true, ...
                    'opticsToEmploy', opticsToEmploy);
        end
    end

    % Stage 3: Optimize pooling from the cone mosaic using the computed
    % cone mosaic STF responses and desired visual RF properties for the
    % RGC mosaic
    if (operationSetToPerformContains.optimizeRGCMosaic)
        % Instantiate the mosaic pooling optimizer
        theMosaicPoolingOptimizer = MosaicPoolingOptimizer(...
            theMidgetRGCMosaic, ...
            'samplingScheme', 'rectangular', ...
            'generateSamplingGrids', true, ...
            'visualizeSamplingGrids', true);

        % Fitting options
        multiStartsNumRetinalPooling = 12;
        multiStartsNumDoGFit = 128;

        % More weight for matching the Rs/Rc ratio
        rmseWeightForRsRcResidual = 2.0;
        rmseWeightForSCintSensResidual = 1.0;

        if (~isempty(arbitraryNodesToCompute))&&(iscell(arbitraryNodesToCompute))
            % Doing an arbitrary selection of nodes (L or M)
            allGridNodesToCompute = arbitraryNodesToCompute;
        else
            % Doing all even, all odd, or all nodes (both L- and M-)
            allGridNodesToCompute = 1:theMosaicPoolingOptimizer.gridNodesNum;
        end

        % Query user which node indices to compute
        nodeIndicesToCompute = input('Which node indices set to compute ? (even or odd) : ', 's');
        switch nodeIndicesToCompute
            case 'even'
                gridNodesToCompute = allGridNodesToCompute(2:2:numel(allGridNodesToCompute));
            case 'odd'
                gridNodesToCompute = allGridNodesToCompute(1:2:numel(allGridNodesToCompute));
            case 'all'
                gridNodesToCompute = allGridNodesToCompute(1:1:numel(allGridNodesToCompute));
            otherwise
                error('grid nodes must be either ''even'' or ''odd''.');
        end
        

        % Compute optimized RGC models (surround cone pooling weights)
        % for each grid node in the RGC mosaic
        for iNode = 1:numel(gridNodesToCompute)

            if (iscell(gridNodesToCompute))
                gridNodeIndex = gridNodesToCompute{iNode}.number;
                whichConeType = gridNodesToCompute{iNode}.coneType;
            else
                gridNodeIndex = gridNodesToCompute(iNode);
                whichConeType = [cMosaic.LCONE_ID cMosaic.MCONE_ID];
            end
            
           
            theMosaicPoolingOptimizer.compute(gridNodeIndex, whichConeType, opticsParams, ...
                fullfile(resourcesDirectory, coneMosaicSTFresponsesFileName), ...
                fullfile(resourcesDirectory, optimizedRGCpoolingObjectsFileName), ...
                'multiStartsNumDoGFit', multiStartsNumDoGFit, ...
                'multiStartsNumRetinalPooling', multiStartsNumRetinalPooling, ...
                'rmseWeightForRsRcResidual', rmseWeightForRsRcResidual, ...
                'rmseWeightForSCintSensResidual', rmseWeightForSCintSensResidual, ...
                'retinalRFmodelParams', retinalRFmodelParams, ...
                'displayFittingProgress', true);

        end
    end

    % Stage 4: Inspect optimized RGC models
    if (operationSetToPerformContains.inspectOptimizedRGCmodels)
        % Instantiate the mosaic pooling optimizer
        theMosaicPoolingOptimizer = MosaicPoolingOptimizer(...
            theMidgetRGCMosaic, ...
            'samplingScheme', 'rectangular', ...
            'generateSamplingGrids', true, ...
            'visualizeSamplingGrids', true);

        gridNodesToInspect = input('Enter grid node to inspect. Hit enter to inspect all.: ');
        if (isempty(gridNodesToInspect))
            gridNodesToInspect = 1:theMosaicPoolingOptimizer.gridNodesNum;
        end

        for iNode = 1:numel(gridNodesToInspect)
            gridNodeIndex = gridNodesToInspect(iNode);
            theMosaicPoolingOptimizer.inspect(gridNodeIndex, opticsParams, ...
                fullfile(resourcesDirectory, optimizedRGCpoolingObjectsFileName));
        end
    end


%     if (operationSetToPerformContains.animateModelConvergence)
%         % Instantiate the mosaic pooling optimizer
%         theMosaicPoolingOptimizer = MosaicPoolingOptimizer(...
%             theMidgetRGCMosaic, ...
%             'generateSamplingGrids', true);
% 
%         gridNodeToAnimate = input('Enter grid node to animate: ');
%         theMosaicPoolingOptimizer.animateModelConvergence(gridNodeToAnimate, ...
%                 fullfile(resourcesDirectory, optimizedRGCpoolingObjectsFileName));
%     end


    % Stage 5: Bake in the computed RGC models (across all grid nodes)
    % to the midgetRGCmosaic, to generate a compute-ready 
    if (operationSetToPerformContains.generateComputeReadyMidgetRGCMosaic)
        % Instantiate the mosaic pooling optimizer with theMidgetRGCMosaic
        theMosaicPoolingOptimizer = MosaicPoolingOptimizer(...
            theMidgetRGCMosaic, ...
            'samplingScheme', 'rectangular', ...
            'generateSamplingGrids', true, ...
            'visualizeSamplingGrids', true);

        visualizeInterpolation = input('Visualize interpolation? [y = YES] : ', 's');
        if (strcmp(visualizeInterpolation, 'y'))
            visualizeInterpolation = true;
        else
            visualizeInterpolation = false;
        end
        
        theMosaicPoolingOptimizer.generateComputeReadyMidgetRGCMosaic(...
            fullfile(resourcesDirectory, optimizedRGCpoolingObjectsFileName), ...
            fullfile(resourcesDirectory, computeReadyMosaicFileName), ...
            visualizeInterpolation);
    end


    % Stage 6: Compute the visual STFs of all cells in the generated compute-ready mRGCMosaic
    if (operationSetToPerformContains.computeVisualSTFsAcrossTheComputeReadyMidgetRGCMosaic)
        MosaicPoolingOptimizer.computeVisualSTFsOfComputeReadyMidgetRGCMosaic(...
            fullfile(resourcesDirectory, computeReadyMosaicFileName), ...
            fullfile(resourcesDirectory, coneMosaicSTFresponsesFileName), ...
            fullfile(resourcesDirectory, mRGCMosaicSTFresponsesFileName));
    end


    % Stage 7: Visualize spatialRFs of the compute-ready mosaic
    if (operationSetToPerformContains.visualizeSpatialRFsAcrossTheComputeReadyMidgetRGCMosaic)
        
        targetRGCposition = input('Enter (xy) position of target RGC (e.g., [5.6 -1.3]): ');
        targetCenterConesNum = input('Enter # of center cones num (e.g, 3): ');
        targetCenterConeMajorityType = input('Enter type of majority center cone num (either cMosaic.LCONE_ID or cMosaic.MCONE_ID): ');

        [~,~,pdfDirectory] = MosaicPoolingOptimizer.resourceFileNameAndPath('pdfsDirectory');
        MosaicPoolingOptimizer.visualizeSpatialRFsAcrossTheComputeReadyMidgetRGCMosaic(...
            fullfile(resourcesDirectory, computeReadyMosaicFileName), ...
            fullfile(resourcesDirectory, mRGCMosaicSTFresponsesFileName), ...
            fullfile(pdfDirectory, 'spatialRFmaps.pdf'), ...
            targetRGCposition, targetCenterConesNum, targetCenterConeMajorityType);
    end

    % Stage 8: Fit the DoG model to the computed visual STFs of all cells in the generated compute-ready mRGCMosaic
    if (operationSetToPerformContains.fitVisualSTFsAcrossTheComputeReadyMidgetRGCMosaic)
        MosaicPoolingOptimizer.fitVisualSTFsOfComputeReadyMidgetRGCMosaic(...
            fullfile(resourcesDirectory, computeReadyMosaicFileName), ...
            fullfile(resourcesDirectory, mRGCMosaicSTFresponsesFileName));
    end

    % Stage 9: Visualized the fitted DoG models to the computed visual STFs of all cells in the generated compute-ready mRGCMosaic
    if (operationSetToPerformContains.visualizeFittedSTFsAcrossTheComputeReadyMidgetRGCMosaic)

        theMosaicFileNames = {...
            fullfile(resourcesDirectory, computeReadyMosaicFileName) ...
            };
        theMRGCSTFResponsesFileNames = { ...
            fullfile(resourcesDirectory, mRGCMosaicSTFresponsesFileName) ...
            };

        employTemporalEquivalentEccentricity = ~true;
        MosaicPoolingOptimizer.visualizeFittedSTFsOfComputeReadyMidgetRGCMosaic(...
            theMosaicFileNames, ...
            theMRGCSTFResponsesFileNames, ...
            employTemporalEquivalentEccentricity);

    end

    

    % Stage 10: Visualized the fitted DoG models to the computed visual STFs of all cells in multiple compute-ready mRGCMosaics
    if (operationSetToPerformContains.visualizeFittedSTFsAcrossMultipleComputeReadyMidgetRGCMosaics)

        mosaicEccs = [2.5 7.0];
        theMosaicFileNames = cell(1, numel(mosaicEccs));
        theMRGCSTFResponsesFileNames = cell(1, numel(mosaicEccs));
        employTemporalEquivalentEccentricity = ~true;

        for iMosaic = 1:numel(mosaicEccs)
            mosaicParams = getMosaicParams(mosaicEccs(iMosaic));

            % Generate the filename of the compute-ready mRGCMosaic
            [computeReadyMosaicFileName, resourcesDirectory] = ...
            MosaicPoolingOptimizer.resourceFileNameAndPath('computeReadyMosaic', ...
                'mosaicParams', mosaicParams, ...
                'opticsParams', opticsParams, ...
                'retinalRFmodelParams', retinalRFmodelParams);

            % Generate the filename of the  mRGCMosaic STF responses
            [mRGCMosaicSTFresponsesFileName, resourcesDirectory] = ...
            MosaicPoolingOptimizer.resourceFileNameAndPath('mRGCMosaicSTFresponses', ...
            'mosaicParams', mosaicParams, ...
            'opticsParams', opticsParams);

            theMosaicFileNames{iMosaic} = fullfile(resourcesDirectory, computeReadyMosaicFileName);
            theMRGCSTFResponsesFileNames{iMosaic} = fullfile(resourcesDirectory, mRGCMosaicSTFresponsesFileName);
        end

        MosaicPoolingOptimizer.visualizeFittedSTFsOfComputeReadyMidgetRGCMosaic(...
            theMosaicFileNames, ...
            theMRGCSTFResponsesFileNames, ...
            employTemporalEquivalentEccentricity);
    end


    % Stage 11: Compute visual RFs
    if (operationSetToPerformContains.computeVisualRFsAcrossTheComputeReadyMidgetRGCMosaic)

        % Load the compute-ready MRGC mosaic
        load(fullfile(resourcesDirectory, computeReadyMosaicFileName), 'theComputeReadyMRGCmosaic');

        % Generate native optics
        theComputeReadyMRGCmosaic.generateNativeOptics(opticsParams);

        reComputeInputConeMosaicSubspaceRFmappingResponses = ~true
        reComputeMRGCMosaicSubspaceRFmappingResponses = true;
        fprintf(2, '\nHit enter to continue:')
        pause
        MosaicPoolingOptimizer.computeVisualRFsOfComputeReadyMidgetRGCMosaic(...
            theComputeReadyMRGCmosaic, ...
            fullfile(resourcesDirectory, coneMosaicSubspaceResponsesFileName), ...
            fullfile(resourcesDirectory, mRGCMosaicSubspaceRresponsesFileName), ...
            reComputeInputConeMosaicSubspaceRFmappingResponses, ...
            reComputeMRGCMosaicSubspaceRFmappingResponses)
    end




%         theROI = regionOfInterest('shape', 'line', 'from', [1 0], 'to', [2.5 0], 'thickness', 0.01);
%         samplingPoints = 1000;  % sample the perimeter of the ROI along 1000 points');
%         pointsPerSample = 4;  % return up to 6 points for each sample along the perimeter');
%         maxDistance = 0.01;     % points must be no further than 0.03 degs away from the closest perimeter sample');
% 
%         % Find cones around the ROI
%         idx = theROI.indicesOfPointsAround(...
%              theMidgetRGCMosaic.inputConeMosaic.coneRFpositionsDegs, pointsPerSample, samplingPoints, maxDistance);
%         identifiedConeIndices = idx;
%         numel(idx)
%         load(inputConeMosaicSTFresponsesFileName);
%         n = reshape(theConeMosaicNullResponses, [1 1 1 numel(theConeMosaicNullResponses)]);
%         
%         for iSF = 1:size(theConeMosaicSTFresponses,2)
%             
%             r = theConeMosaicSTFresponses(7,iSF,:,:);
%             
%             c = squeeze(bsxfun(@times, bsxfun(@minus,r,n),1./n));
% size(c)
%     figure()
%     c =c(:,idx);
%     [min(c(:)) max(c(:))]
%             STF = max(c,[],1);
%             imagesc(c, [-1 1])
% colormap(gray(1024))
%             max(STF(:))
%             size(STF)
%             figure(99); clf;
%             lConeIndices = find(theMidgetRGCMosaic.inputConeMosaic.coneTypes(idx) == cMosaic.LCONE_ID)
%             mConeIndices = find(theMidgetRGCMosaic.inputConeMosaic.coneTypes(idx) == cMosaic.MCONE_ID)
%             plot(theMidgetRGCMosaic.inputConeMosaic.coneRFpositionsDegs(idx(lConeIndices),1), STF(lConeIndices), 'rs-')
%             hold on;
%             plot(theMidgetRGCMosaic.inputConeMosaic.coneRFpositionsDegs(idx(mConeIndices),1), STF(mConeIndices), 'gs-')
%             pause
%             drawnow;
%         end

end


function mosaicParams = getMosaicParams(mosaicEcc)

    switch (mosaicEcc)
        case 0
            % Mosaic params to employ. This is for the 2.5 deg - centered mosaic
        % which covers the [1 - 4] deg eccentricity range
        mosaicParams = struct(...
            'eccDegs', [0 0], ...
            'sizeDegs', [3 3]);

        case 2.5
        % Mosaic params to employ. This is for the 2.5 deg - centered mosaic
        % which covers the [1 - 4] deg eccentricity range
        mosaicParams = struct(...
            'eccDegs', [2.5 0], ...
            'sizeDegs', [3 3]);

        case 7.0
        % Mosaic params to employ. This is for the 7.0 deg - centered mosaic
        % which covers the [4-10] deg eccentricity range
        mosaicParams = struct(...
            'eccDegs', [7 0], ...
            'sizeDegs', [6 3]);

        case -10.0
        mosaicParams = struct(...
            'eccDegs', [-10 0], ...
            'sizeDegs', [6 3]);



        otherwise
            error('No data for this eccentricity')
    end
end


function arbitraryNodesToCompute = selectNodesToRecompute()
    arbitraryNodesToCompute = {};

    s = input('Compute all nodes (1) or speficic ones (2):');
    if (s == 1)
        arbitraryNodesToCompute
        return;
    else

        % Bad nodes
        arbitraryNodesToCompute{numel(arbitraryNodesToCompute)+1} = struct(...
                'number', 5, ...
                'coneType', [cMosaic.LCONE_ID cMosaic.MCONE_ID]);
    end
end


function operationSetToPerformContains = promptUserForOperationsToPerform(mosaicParams)

    % Generate the center-connected mosaic
    operationSetToPerformContains.generateRGCMosaic = ~true;

    % Visualize the center-connected mosaic
    operationSetToPerformContains.visualizeRGCMosaic = ~true;

    % Visualize PSFs at various positions within the mosaic
    operationSetToPerformContains.visualizePSFsAtEccentricities = ~true;

    % Compute cone mosaic STF responses
    operationSetToPerformContains.computeConeMosaicSTFresponses = ~true;

    % Derive optimized surround cone pooling kernels (takes a long time)
    operationSetToPerformContains.optimizeRGCMosaic = ~true;

    % Examine optimized cone pooling models at all grids
    operationSetToPerformContains.inspectOptimizedRGCmodels = ~true;

    % Generate the compute-ready mRGC mosaic
    operationSetToPerformContains.generateComputeReadyMidgetRGCMosaic = ~true;

    % Visualize the derived spatial RFs of the compute-ready mRGC mosaic
    operationSetToPerformContains.visualizeSpatialRFsAcrossTheComputeReadyMidgetRGCMosaic = ~true;

    % Validate the compute-ready mRGC mosaic
    operationSetToPerformContains.computeVisualSTFsAcrossTheComputeReadyMidgetRGCMosaic = ~true;
    operationSetToPerformContains.fitVisualSTFsAcrossTheComputeReadyMidgetRGCMosaic = ~true;
    operationSetToPerformContains.visualizeFittedSTFsAcrossTheComputeReadyMidgetRGCMosaic = ~true;
    operationSetToPerformContains.visualizeFittedSTFsAcrossMultipleComputeReadyMidgetRGCMosaics = ~true;
    operationSetToPerformContains.computeVisualRFsAcrossTheComputeReadyMidgetRGCMosaic = ~true;
    

    operationSetToPerformContains.animateModelConvergence = ~true;

    actionStrings{1} = '[ 1] Generate the center-connected mRGCMosaic';
    actionStrings{2} = '[ 2] Visualize the center-connected mRGCMosaic';
    actionStrings{3} = '[ 3] Visualize PSFs at a 3x3 grid within the mRGCMosaic';
    actionStrings{4} = '[ 4] Compute responses across spatial frequencies and orientations of the input cone mosaic';
    actionStrings{5} = '[ 5] Derive optimized surround cone pooling kernels. (NOTE: This step takes a long time';
    actionStrings{6} = '[ 6] Examine optimized cone pooling models at all or certain grid nodes within the mosaic';
    actionStrings{7} = '[ 7] Generate a compute-ready (center-surround connected) mRGCMosaic\n\t   based on the previously derived optimized surround cone pooling kernels';
    actionStrings{8} = '[ 8] Compute-ready mRGCMosaic validation. Step1: compute visual STFs for all cells in the mosaic';
    actionStrings{9} = '[ 9] Compute-ready mRGCMosaic validation. Step2: fit a DoG model to the computed visual STFs for all cells in the mosaic';
    actionStrings{10} = '[10] Compute-ready mRGCMosaic: visualize spatial RFs and visual STF';
    actionStrings{11} = '[11] Compute-ready mRGCMosaic validation. Step3: visualize fitted DoG model params for all cells in the mosaic';
    actionStrings{12} = '[12] Compute-ready mRGCMosaic validation. Step4: visualize fitted DoG model params for all cells in multiple mosaics';
    actionStrings{13} = '[13] Compute-ready mRGCMosaic: compute visual RFs for all cells in the mosaic';

    invalidActionSelected = true;
    validChoiceIDs = 1:numel(actionStrings);
    while (invalidActionSelected)
        fprintf('\n\nAvailable actions for the mosaic at eccentricity (%2.1f,%2.1f):', ...
            mosaicParams.eccDegs(1), mosaicParams.eccDegs(2));
        for iString = 1:numel(actionStrings)
            fprintf('\n\n\t%s', actionStrings{iString});
        end
        choice = input('\n\nEnter action ID : ', 's');
        if (~isempty(choice))
            choiceID = str2num(choice);
            if (ismember(choiceID, validChoiceIDs))
                switch (choiceID)
                    case 1
                        % Generate center-connected mosaic
                        operationSetToPerformContains.generateRGCMosaic = true;
                    case 2
                        % Visualize center-connected mosaic
                        operationSetToPerformContains.visualizeRGCMosaic = true;
                    case 3
                        % Visualize PSFs at various positions within the mosaic
                        operationSetToPerformContains.visualizePSFsAtEccentricities = true;
                    case 4
                        % Compute cone mosaic STF responses
                        operationSetToPerformContains.computeConeMosaicSTFresponses = true;
                    case 5
                        % Derive optimized surround cone pooling kernels (takes a long time)
                        operationSetToPerformContains.optimizeRGCMosaic = true;
                    case 6
                        % Examine optimized cone pooling models at all or
                        % certain  grid nodes
                        operationSetToPerformContains.inspectOptimizedRGCmodels = true;
                    case 7
                        % Generate a compute-ready MRGC mosaic
                        operationSetToPerformContains.generateComputeReadyMidgetRGCMosaic = true;
                    case 8
                        % Validate a compute-ready mRGCMosaic: step1 - compute visual STFs for all cells
                        operationSetToPerformContains.computeVisualSTFsAcrossTheComputeReadyMidgetRGCMosaic = true;
                    case 9
                        % Validate a compute-ready mRGCMosaic: step2 - fit a DoG model to all the computed visual STFs for all cells
                        operationSetToPerformContains.fitVisualSTFsAcrossTheComputeReadyMidgetRGCMosaic = true;
                    case 10
                        % Visualize the derived spatial RFs of the compute-ready mRGC mosaic
                        operationSetToPerformContains.visualizeSpatialRFsAcrossTheComputeReadyMidgetRGCMosaic = true;
                    case 11
                        % Validate a compute-ready mRGCMosaic: step3 - visualize fitted visual STFs for all cells in a single mRGC mosaic
                        operationSetToPerformContains.visualizeFittedSTFsAcrossTheComputeReadyMidgetRGCMosaic = true;
                    case 12
                        % Validate a compute-ready mRGCMosaic: step4 - visualize fitted visual STFs for all cells in multiple mRGC mosaics
                        operationSetToPerformContains.visualizeFittedSTFsAcrossMultipleComputeReadyMidgetRGCMosaics = true;
                    case 13
                        operationSetToPerformContains.computeVisualRFsAcrossTheComputeReadyMidgetRGCMosaic = true;
                    otherwise
                        error('Unknown option')

                end % switch
                invalidActionSelected = false;
            end
        end

    end
   
end
