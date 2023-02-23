function testMidgetRGCMosaic()

    
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
    optimizeRGCMosaic = true;

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


    opticsParams = theMidgetRGCMosaic.defaultOpticsParams;
    [coneMosaicSTFresponsesFileName, resourcesDirectory] = ...
        MosaicPoolingOptimizer.resourceFileNameAndPath('coneMosaicSTFresponses', ...
            'mosaicParams', mosaicParams, ...
            'opticsParams', opticsParams);
    
    [optimizedRGCpoolingObjectsFileName, resourcesDirectory] = ...
        MosaicPoolingOptimizer.resourceFileNameAndPath('optimizedRGCpoolingObjects', ...
            'mosaicParams', mosaicParams, ...
            'opticsParams', opticsParams);

    
    % Stage 2: Generate optics and compute input cone mosaic STF responses
    if (computeConeMosaicSTFresponses)

        % Generate the native optics. These will be used to compute the
        % cone mosaic STF responses and optimize the RGC RFs
        theMidgetRGCMosaic.generateNativeOptics(opticsParams);

        

        compuseConeMosaicSTFResponsesSeparatelyAtEachGridLocation = false;
        if (compuseConeMosaicSTFResponsesSeparatelyAtEachGridLocation)
            % Compute the input cone mosaic visual STFs, at a single
            % gridNode (using small stimulus patches centered on that node)

            % Instantiate the mosaic pooling optimizer
            theMosaicPoolingOptimizer = MosaicPoolingOptimizer(...
                theMidgetRGCMosaic, ...
                'generateSamplingGrids', true);

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
        % Instantiate the mosaic pooling optimizer
        theMosaicPoolingOptimizer = MosaicPoolingOptimizer(...
            theMidgetRGCMosaic, ...
            'generateSamplingGrids', true);

        % Compute optimized RGC models (surround cone pooling weights)
        % for each grid node in the RGC mosaic
        for gridNodeIndex = 1:theMosaicPoolingOptimizer.gridNodesNum

            theMosaicPoolingOptimizer.compute(gridNodeIndex, ...
                fullfile(resourcesDirectory, coneMosaicSTFresponsesFileName), ...
                fullfile(resourcesDirectory, optimizedRGCpoolingObjectsFileName), ...
                'multiStartsNumDoGFit', 64, ...
                'multiStartsNumRetinalPooling', 8, ...
                'displayFittingProgress', true);

        end

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
