function testMidgetRGCMosaic
    
    arbitraryNodesToCompute = []; % selectNodesToRecompute();

    % Mosaic params to employ
    mosaicParams = struct(...
        'eccDegs', [2.5 0], ...
        'sizeDegs', [3 3]);

    % Get mosaic filename
    [mosaicFileName, resourcesDirectory] = ...
        MosaicPoolingOptimizer.resourceFileNameAndPath('mosaic', ...
            'mosaicParams', mosaicParams);
    
    % Actions
    generateRGCMosaic = ~true;
    computeConeMosaicSTFresponses = ~true;
    optimizeRGCMosaic = ~true;
    inspectOptimizedRGCmodels = ~true;
    generateComputeReadyMidgetRGCMosaic = ~true;
    computeVisualSTFsAcrossTheComputeReadyMidgetRGCMosaic =~true;
    inspectVisualSTFsAcrossTheComputeReadyMidgetRGCMosaic = true;
    animateModelConvergence = ~true;

    if (generateRGCMosaic)
        % Generate mosaic, its input coneMosaic and connect cones to the RF centers
        theMidgetRGCMosaic = mRGCMosaic(...
            'eccentricityDegs', mosaicParams.eccDegs, ...
            'sizeDegs', mosaicParams.sizeDegs);
    
        save(fullfile(resourcesDirectory,mosaicFileName), 'theMidgetRGCMosaic');
    else
        load(fullfile(resourcesDirectory,mosaicFileName), 'theMidgetRGCMosaic');
    end
    
    visualizeMosaics = false;
    if (visualizeMosaics)
        % Visualize input cone mosaic
        theMidgetRGCMosaic.inputConeMosaic.visualize()
    
        % Visualize cone pooling by the RF centers
        theMidgetRGCMosaic.visualize();
    end

   
    % Optics params to employ
    opticsParams = theMidgetRGCMosaic.defaultOpticsParams;

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

    % Stage 2: Generate optics and compute input cone mosaic STF responses
    if (computeConeMosaicSTFresponses)

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
    if (optimizeRGCMosaic)
        % Make sure 
        nodeIndicesToCompute = input('Which node indices set to compute ? (even or odd) : ', 's');

        
        assert(...
            (ismember(nodeIndicesToCompute, {'even', 'odd', 'all'})), ...
            'grid nodes indices must be either ''even'' or ''odd''.');

        % Instantiate the mosaic pooling optimizer
        theMosaicPoolingOptimizer = MosaicPoolingOptimizer(...
            theMidgetRGCMosaic, ...
            'generateSamplingGrids', true);

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
        

        % Generate filename for the optimizedRGCpoolingObjects
        [optimizedRGCpoolingObjectsFileName, resourcesDirectory] = ...
        MosaicPoolingOptimizer.resourceFileNameAndPath('optimizedRGCpoolingObjects', ...
            'mosaicParams', mosaicParams, ...
            'opticsParams', opticsParams, ...
            'retinalRFmodelParams', retinalRFmodelParams);


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
    if (inspectOptimizedRGCmodels)
        % Instantiate the mosaic pooling optimizer
        theMosaicPoolingOptimizer = MosaicPoolingOptimizer(...
            theMidgetRGCMosaic, ...
            'generateSamplingGrids', true);

        gridNodesToInspect = input('Enter grid node to inspect. Hit enter to inspect all.: ');
        if (isempty(gridNodesToInspect))
            gridNodesToInspect = 1:theMosaicPoolingOptimizer.gridNodesNum;
        end
        

        [optimizedRGCpoolingObjectsFileName, resourcesDirectory] = ...
        MosaicPoolingOptimizer.resourceFileNameAndPath('optimizedRGCpoolingObjects', ...
            'mosaicParams', mosaicParams, ...
            'opticsParams', opticsParams, ...
            'retinalRFmodelParams', retinalRFmodelParams);

        for iNode = 1:numel(gridNodesToInspect)
            gridNodeIndex = gridNodesToInspect(iNode);
            theMosaicPoolingOptimizer.inspect(gridNodeIndex, opticsParams, ...
                fullfile(resourcesDirectory, optimizedRGCpoolingObjectsFileName));
        end

    end

    if (animateModelConvergence)
        % Instantiate the mosaic pooling optimizer
        theMosaicPoolingOptimizer = MosaicPoolingOptimizer(...
            theMidgetRGCMosaic, ...
            'generateSamplingGrids', true);

        [optimizedRGCpoolingObjectsFileName, resourcesDirectory] = ...
        MosaicPoolingOptimizer.resourceFileNameAndPath('optimizedRGCpoolingObjects', ...
            'mosaicParams', mosaicParams, ...
            'opticsParams', opticsParams, ...
            'retinalRFmodelParams', retinalRFmodelParams);

        gridNodeToAnimate = input('Enter grid node to animate: ');
        theMosaicPoolingOptimizer.animateModelConvergence(gridNodeToAnimate, ...
                fullfile(resourcesDirectory, optimizedRGCpoolingObjectsFileName));
    end


    % Stage 5: Bake in the computed RGC models (across all grid nodes)
    % to the midgetRGCmosaic, to generate a compute-ready 
    if (generateComputeReadyMidgetRGCMosaic)
        % Instantiate the mosaic pooling optimizer with theMidgetRGCMosaic
        theMosaicPoolingOptimizer = MosaicPoolingOptimizer(...
            theMidgetRGCMosaic, ...
            'generateSamplingGrids', true);

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

        visualizeInterpolation = input('Visualize interpolation? [y = YES] : ', 's');
        if (strcmp(visualizeInterpolation, 'y'))
            visualizeInterpolation = true;
        else
            visualizeInterpolation = false;
        end
        
        theMosaicPoolingOptimizer.generateComputeReadyMidgetRGCMosaic(...
            fullfile(resourcesDirectory, optimizedRGCpoolingObjectsFileName), ...
            fullfile(resourcesDirectory, computeReadyMosaicFileName), ...
            visualizeInterpolation)

    end


    % Stage 6: Compute the visual STFs of all cells in the generated compute-ready mRGCMosaic
    if (computeVisualSTFsAcrossTheComputeReadyMidgetRGCMosaic)

        % Generate the filename of the compute-ready mRGCMosaic to generate
        [computeReadyMosaicFileName, resourcesDirectory] = ...
            MosaicPoolingOptimizer.resourceFileNameAndPath('computeReadyMosaic', ...
                'mosaicParams', mosaicParams, ...
                'opticsParams', opticsParams, ...
                'retinalRFmodelParams', retinalRFmodelParams);

        MosaicPoolingOptimizer.computeVisualSTFsOfComputeReadyMidgetRGCMosaic(...
            fullfile(resourcesDirectory, computeReadyMosaicFileName), ...
            fullfile(resourcesDirectory, coneMosaicSTFresponsesFileName), ...
            fullfile(resourcesDirectory, mRGCMosaicSTFresponsesFileName));
    end


    % Stage 6: Fit the DoG model to the computed visual STFs of all cells in the generated compute-ready mRGCMosaic
    if ( inspectVisualSTFsAcrossTheComputeReadyMidgetRGCMosaic)

        % Generate the filename of the compute-ready mRGCMosaic to generate
        [computeReadyMosaicFileName, resourcesDirectory] = ...
            MosaicPoolingOptimizer.resourceFileNameAndPath('computeReadyMosaic', ...
                'mosaicParams', mosaicParams, ...
                'opticsParams', opticsParams, ...
                'retinalRFmodelParams', retinalRFmodelParams);

        MosaicPoolingOptimizer.inspectVisualSTFsOfComputeReadyMidgetRGCMosaic(...
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


function arbitraryNodesToCompute = selectNodesToRecompute()
    arbitraryNodesToCompute = {};

    arbitraryNodesToCompute{numel(arbitraryNodesToCompute)+1} = struct(...
        'number', 5, ...
        'coneType', [cMosaic.LCONE_ID]);

    arbitraryNodesToCompute{numel(arbitraryNodesToCompute)+1} = struct(...
        'number', 7, ...
        'coneType', [cMosaic.MCONE_ID]);

    arbitraryNodesToCompute{numel(arbitraryNodesToCompute)+1} = struct(...
        'number', 13, ...
        'coneType', [cMosaic.MCONE_ID]);

    arbitraryNodesToCompute{numel(arbitraryNodesToCompute)+1} = struct(...
        'number', 14, ...
        'coneType', [cMosaic.LCONE_ID]);


end
