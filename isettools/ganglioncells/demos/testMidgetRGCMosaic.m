function testMidgetRGCMosaic
    
    arbitraryNodesToCompute =  selectNodesToRecompute();
    mosaicEcc = 2.5;
    mosaicEcc = 7.0;
    
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
        % Generate the native optics appropriate for the midgetRGCMosaic at hand. 
        % These will be used to compute the
        % cone mosaic STF responses and optimize the RGC RFs
        theMidgetRGCMosaic.generateNativeOptics(opticsParams);

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
                theMosaicPoolingOptimizer.generateConeMosaicSTFresponses(...
                    gridNodeIndex, stimSizeDegs,  ...
                    fullfile(resourcesDirectory, coneMosaicSTFresponsesFileNameForThisNode), ...
                    'useParfor', true, ...
                    'visualizedResponses', ~true);
            end
        else
            % Instantiate the mosaic pooling optimizer
            theMosaicPoolingOptimizer = MosaicPoolingOptimizer(...
                theMidgetRGCMosaic, ...
                'generateSamplingGrids', false);

            % Compute the input cone mosaic visual STFs across all nodes, 
            % using full field stimuli
            theMosaicPoolingOptimizer.generateConeMosaicSTFresponses(...
                    [], [], ...
                    fullfile(resourcesDirectory, coneMosaicSTFresponsesFileName), ...
                    'useParfor', ~true, ...
                    'visualizedResponses', ~true);
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


    % Stage 6: Fit the DoG model to the computed visual STFs of all cells in the generated compute-ready mRGCMosaic
    if (operationSetToPerformContains.fitVisualSTFsAcrossTheComputeReadyMidgetRGCMosaic)
        MosaicPoolingOptimizer.fitVisualSTFsOfComputeReadyMidgetRGCMosaic(...
            fullfile(resourcesDirectory, computeReadyMosaicFileName), ...
            fullfile(resourcesDirectory, mRGCMosaicSTFresponsesFileName));
    end

    % Stage 7: Visualized the fitted DoG models to the computed visual STFs of all cells in the generated compute-ready mRGCMosaic
    if (operationSetToPerformContains.visualizeFittedVisualSTFsAcrossTheComputeReadyMidgetRGCMosaic)
        MosaicPoolingOptimizer.visualizeFittedVisualSTFsOfComputeReadyMidgetRGCMosaic(...
            fullfile(resourcesDirectory, computeReadyMosaicFileName), ...
            fullfile(resourcesDirectory, mRGCMosaicSTFresponsesFileName));
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

        for i = 7:27
            arbitraryNodesToCompute{numel(arbitraryNodesToCompute)+1} = struct(...
                    'number', i, ...
                    'coneType', [cMosaic.LCONE_ID cMosaic.MCONE_ID]);
        end
    end
end


function operationSetToPerformContains = promptUserForOperationsToPerform(mosaicParams)

    % Generate mosaic
    operationSetToPerformContains.generateRGCMosaic = ~true;

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

    % Validate the compute-ready mRGC mosaic
    operationSetToPerformContains.computeVisualSTFsAcrossTheComputeReadyMidgetRGCMosaic = ~true;
    operationSetToPerformContains.fitVisualSTFsAcrossTheComputeReadyMidgetRGCMosaic = ~true;
    operationSetToPerformContains.visualizeFittedVisualSTFsAcrossTheComputeReadyMidgetRGCMosaic = ~true;

    operationSetToPerformContains.animateModelConvergence = ~true;

    invalidActionSelected = true;
    validChoiceIDs = 1:9;
    while (invalidActionSelected)
        fprintf('\n\nAvailable actions for the mosaic at eccentricity (%2.1f,%2.1f):', ...
            mosaicParams.eccDegs(1), mosaicParams.eccDegs(2));
        fprintf('\n\n\t[1] Generate the center-connected mRGCMosaic');
        fprintf('\n\n\t[2] Visualize PSFs at a 3x3 grid within the mRGCMosaic');
        fprintf('\n\n\t[3] Compute responses across spatial frequencies and orientations of the input cone mosaic'); 
        fprintf('\n\n\t[4] Derive optimized surround cone pooling kernels. (NOTE: This step takes a long time)');
        fprintf('\n\n\t[5] Examine optimized cone pooling models at all or certain grid nodes within the mosaic');
        fprintf('\n\n\t[6] Generate a compute-ready (center-surround connected) mRGCMosaic\n\t   based on the previously derived optimized surround cone pooling kernels');
        fprintf('\n\n\t[7] Compute-ready mRGCMosaic validation. Step1: compute visual STFs for all cells in the mosaic');
        fprintf('\n\n\t[8] Compute-ready mRGCMosaic validation. Step2: fit a DoG model to the computed visual STFs for all cells in the mosaic');
        fprintf('\n\n\t[9] Compute-ready mRGCMosaic validation. Step3: visualize fitted DoG model params for all cells in the mosaic');

        choice = input('\n\nEnter action ID : ', 's');
        if (~isempty(choice))
            choiceID = str2num(choice);
            if (ismember(choiceID, validChoiceIDs))
                switch (choiceID)
                    case 1
                        % Generate mosaic
                        operationSetToPerformContains.generateRGCMosaic = true;
                    case 2
                        % Visualize PSFs at various positions within the mosaic
                        operationSetToPerformContains.visualizePSFsAtEccentricities = true;
                    case 3
                        % Compute cone mosaic STF responses
                        operationSetToPerformContains.computeConeMosaicSTFresponses = true;
                    case 4
                        % Derive optimized surround cone pooling kernels (takes a long time)
                        operationSetToPerformContains.optimizeRGCMosaic = true;
                    case 5
                        % Examine optimized cone pooling models at all or
                        % certain  grid nodes
                        operationSetToPerformContains.inspectOptimizedRGCmodels = true;
                    case 6
                        % Generate a compute-ready MRGC mosaic
                        operationSetToPerformContains.generateComputeReadyMidgetRGCMosaic = true;
                    case 7
                        % Validate a compute-ready mRGCMosaic: step1 - compute visual STFs for all cells
                        operationSetToPerformContains.computeVisualSTFsAcrossTheComputeReadyMidgetRGCMosaic = true;
                    case 8
                        % Validate a compute-ready mRGCMosaic: step2 - fit a DoG model to all the computed visual STFs for all cells
                        operationSetToPerformContains.fitVisualSTFsAcrossTheComputeReadyMidgetRGCMosaic = true;
                    case 9
                        % Validate a compute-ready mRGCMosaic: step3 - visualize fitted visual STFs for all cells
                        operationSetToPerformContains.visualizeFittedVisualSTFsAcrossTheComputeReadyMidgetRGCMosaic = true;
                    otherwise
                        error('Unknown option')

                end % switch
                invalidActionSelected = false;
            end
        end

    end
   
end
