function testMidgetRGCMosaic

    
    arbitraryNodesToCompute = [] %selectNodesToRecompute();

%     arbitraryNodesToCompute{1} = struct(...
%         'number', 16, ...
%         'coneType', [cMosaic.LCONE_ID]);

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
    visualizeConeMosaicSTFresponses = ~true;
    inspectOptimizedRGCmodels = true;
    generateComputeReadyMidgetRGCMosaic = ~true;


    if (generateRGCMosaic)
        % Generate mosaic, its input coneMosaic and connect cones to the RF centers
        theMidgetRGCMosaic = mRGCMosaic(...
            'eccentricityDegs', mosaicParams.eccDegs, ...
            'sizeDegs', mosaicParams.sizeDegs);
    
        save(fullfile(resourcesDirectory,mosaicFileName), 'theMidgetRGCMosaic');
    else
        load(fullfile(resourcesDirectory,mosaicFileName), 'theMidgetRGCMosaic');
    end 
    
    % Optics params to employ
    opticsParams = theMidgetRGCMosaic.defaultOpticsParams;


    visualizeMosaics = ~true;
    if (visualizeMosaics)

        % Visualize input cone mosaic
        theMidgetRGCMosaic.inputConeMosaic.visualize()
    
        % Visualize cone pooling by the RF centers
        theMidgetRGCMosaic.visualize();
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

    if (visualizeConeMosaicSTFresponses)
        targetConePositions = [...
            0.77 0.0; ...
            1.4 0.0; ...
            2.42 0.0; ...
            3.3 0.0; ...
            4.2 0.0];

        MosaicPoolingOptimizer.visualizeConeMosaicSTFresponses(...
                    fullfile(resourcesDirectory,mosaicFileName), ...
                    fullfile(resourcesDirectory, coneMosaicSTFresponsesFileName), ...
                    'opticsParams', opticsParams, ...
                    'targetConePositions', targetConePositions);
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
            
            theMosaicPoolingOptimizer.compute(gridNodeIndex, whichConeType, ...
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
            'generateSamplingGrids', true, ...
            'visualizeSamplingGrids', true);

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
            theMosaicPoolingOptimizer.inspect(gridNodeIndex, ...
                fullfile(resourcesDirectory, optimizedRGCpoolingObjectsFileName));
        end

    end

    % Stage 5: Bake in the computed RGC models (across all grid nodes)
    % to the midgetRGCmosaic
    if (generateComputeReadyMidgetRGCMosaic)
        % Instantiate the mosaic pooling optimizer with theMidgetRGCMosaic
        theMosaicPoolingOptimizer = MosaicPoolingOptimizer(...
            theMidgetRGCMosaic, ...
            'generateSamplingGrids', true);

        [optimizedRGCpoolingObjectsFileName, resourcesDirectory] = ...
            MosaicPoolingOptimizer.resourceFileNameAndPath('optimizedRGCpoolingObjects', ...
                'mosaicParams', mosaicParams, ...
                'opticsParams', opticsParams, ...
                'retinalRFmodelParams', retinalRFmodelParams);

        [computeReadyMosaicFileName, resourcesDirectory] = ...
            MosaicPoolingOptimizer.resourceFileNameAndPath('computeReadyMosaic', ...
                'mosaicParams', mosaicParams, ...
                'opticsParams', opticsParams, ...
                'retinalRFmodelParams', retinalRFmodelParams);

        theMosaicPoolingOptimizer.generateComputeReadyMidgetRGCMosaic(...
            fullfile(resourcesDirectory, optimizedRGCpoolingObjectsFileName), ...
            fullfile(resourcesDirectory, computeReadyMosaicFileName));

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
        'number', 1, ...
        'coneType', [cMosaic.LCONE_ID cMosaic.MCONE_ID]);

    arbitraryNodesToCompute{numel(arbitraryNodesToCompute)+1} = struct(...
        'number', 3, ...
        'coneType', [cMosaic.MCONE_ID]);

    arbitraryNodesToCompute{numel(arbitraryNodesToCompute)+1} = struct(...
        'number', 4, ...
        'coneType', [cMosaic.LCONE_ID cMosaic.MCONE_ID]);

    arbitraryNodesToCompute{numel(arbitraryNodesToCompute)+1} = struct(...
        'number', 5, ...
        'coneType', [cMosaic.LCONE_ID]);

    arbitraryNodesToCompute{numel(arbitraryNodesToCompute)+1} = struct(...
        'number', 6, ...
        'coneType', [cMosaic.LCONE_ID cMosaic.MCONE_ID]);

    arbitraryNodesToCompute{numel(arbitraryNodesToCompute)+1} = struct(...
        'number', 7, ...
        'coneType', [cMosaic.LCONE_ID cMosaic.MCONE_ID]);

     arbitraryNodesToCompute{numel(arbitraryNodesToCompute)+1} = struct(...
        'number', 9, ...
        'coneType', [cMosaic.LCONE_ID]);

     arbitraryNodesToCompute{numel(arbitraryNodesToCompute)+1} = struct(...
        'number', 10, ...
        'coneType', [cMosaic.LCONE_ID cMosaic.MCONE_ID]);

     arbitraryNodesToCompute{numel(arbitraryNodesToCompute)+1} = struct(...
        'number', 11, ...
        'coneType', [cMosaic.LCONE_ID cMosaic.MCONE_ID]);

     arbitraryNodesToCompute{numel(arbitraryNodesToCompute)+1} = struct(...
        'number', 12, ...
        'coneType', [cMosaic.LCONE_ID]);

     arbitraryNodesToCompute{numel(arbitraryNodesToCompute)+1} = struct(...
        'number', 13, ...
        'coneType', [cMosaic.LCONE_ID cMosaic.MCONE_ID]);

     arbitraryNodesToCompute{numel(arbitraryNodesToCompute)+1} = struct(...
        'number', 15, ...
        'coneType', [cMosaic.LCONE_ID cMosaic.MCONE_ID]);

     arbitraryNodesToCompute{numel(arbitraryNodesToCompute)+1} = struct(...
        'number', 16, ...
        'coneType', [cMosaic.LCONE_ID cMosaic.MCONE_ID]);

     arbitraryNodesToCompute{numel(arbitraryNodesToCompute)+1} = struct(...
        'number', 17, ...
        'coneType', [cMosaic.LCONE_ID cMosaic.MCONE_ID]);

     arbitraryNodesToCompute{numel(arbitraryNodesToCompute)+1} = struct(...
        'number', 18, ...
        'coneType', [cMosaic.LCONE_ID cMosaic.MCONE_ID]);

     arbitraryNodesToCompute{numel(arbitraryNodesToCompute)+1} = struct(...
        'number', 19, ...
        'coneType', [cMosaic.LCONE_ID cMosaic.MCONE_ID]);

     arbitraryNodesToCompute{numel(arbitraryNodesToCompute)+1} = struct(...
        'number', 21, ...
        'coneType', [cMosaic.LCONE_ID]);

     arbitraryNodesToCompute{numel(arbitraryNodesToCompute)+1} = struct(...
        'number', 22, ...
        'coneType', [cMosaic.MCONE_ID]);

     arbitraryNodesToCompute{numel(arbitraryNodesToCompute)+1} = struct(...
        'number', 23, ...
        'coneType', [cMosaic.LCONE_ID]);

     arbitraryNodesToCompute{numel(arbitraryNodesToCompute)+1} = struct(...
        'number', 24, ...
        'coneType', [cMosaic.MCONE_ID]);

     arbitraryNodesToCompute{numel(arbitraryNodesToCompute)+1} = struct(...
        'number', 25, ...
        'coneType', [cMosaic.LCONE_ID cMosaic.MCONE_ID]);

end
