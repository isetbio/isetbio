function testMidgetRGCMosaic()

    
    mosaicEccDegs = [3 0];
    mosaicSizeDegs = [2 2];

    resourcesDirectory = fullfile(...
                MosaicPoolingOptimizer.localDropboxPath, ...
                'productionMidgetRGCMosaics/MosaicOptimizerResources');

    mosaiFileName = sprintf('mRGCMosaicEcDegs(%2.1f_%2.1f)_SizeDegs(%2.1f_%2.1f).mat', ...
        mosaicEccDegs(1), mosaicEccDegs(2), mosaicSizeDegs(1), mosaicSizeDegs(2));

    % Actions
    recomputeMosaic = true;
    recomputeConeMosaicSTFresponses = true;

    if (recomputeMosaic)
        % Generate mosaic, its input coneMosaic and connect cones to the RF centers
        theMidgetRGCMosaic = mRGCMosaic(...
            'eccentricityDegs', mosaicEccDegs, ...
            'sizeDegs', mosaicSizeDegs);
    
        save(fullfile(resourcesDirectory,mosaiFileName), 'theMidgetRGCMosaic');
    else
        load(fullfile(resourcesDirectory,mosaiFileName), 'theMidgetRGCMosaic');
    end

    % Visualize input cone mosaic
    theMidgetRGCMosaic.inputConeMosaic.visualize()

    % Visualize cone pooling by the RF centers
    theMidgetRGCMosaic.visualize();

    % Stage 2: Generate optics and connect cones to the RF surrounds
    if (recomputeConeMosaicSTFresponses)
        % Generate the native optics. These optics will be used 
        % for optimizing surround weights so as to achieve target visual RF
        % properties.
        % Retrieve the default optics params
        opticsParams = theMidgetRGCMosaic.defaultOpticsParams;
        % Optionally, modify certain opticsParams

        % Generate the native optics
        theMidgetRGCMosaic.generateNativeOptics(opticsParams);

        % Instantiate the mosaic pooling optimizer
        theMosaicPoolingOptimizer = MosaicPoolingOptimizer(theMidgetRGCMosaic);

        % Compute cone mosaic STF responses for 1x1 stim patches
        % positioned at the targetRGCposition of each grid node
        for gridNode = 1:theMosaicPoolingOptimizer.gridNodesNum
            stimSizeDegs = [1 1];

            % Responses filename
            responsesFileName = sprintf('STFresponsesGridNode_%d.mat', gridNode);
            responsesFileName = fullfile(resourcesDirectory, responsesFileName );

            theMosaicPoolingOptimizer.computeConeMosaicSTFresponses(...
                gridNode, stimSizeDegs, responsesFileName, ...
                'useParfor', true, ...
                'visualizedResponses', ~true);
        end



        % Fit the cone pooling weight optimizer at a particular grid node
%         theMosaicPoolingOptimizer.fit(...
%             'gridNodeList', [1], ...
%             'recomputeConeMosaicSTFs', true);

    end

    recomputeStage2b = ~true;
    if (recomputeStage2b)

        theROI = regionOfInterest('shape', 'line', 'from', [1 0], 'to', [2.5 0], 'thickness', 0.01);
        samplingPoints = 1000;  % sample the perimeter of the ROI along 1000 points');
        pointsPerSample = 4;  % return up to 6 points for each sample along the perimeter');
        maxDistance = 0.01;     % points must be no further than 0.03 degs away from the closest perimeter sample');

        % Find cones around the ROI
        idx = theROI.indicesOfPointsAround(...
             theMidgetRGCMosaic.inputConeMosaic.coneRFpositionsDegs, pointsPerSample, samplingPoints, maxDistance);
        identifiedConeIndices = idx;
        numel(idx)
        load(inputConeMosaicSTFresponsesFileName);
        n = reshape(theConeMosaicNullResponses, [1 1 1 numel(theConeMosaicNullResponses)]);
        
        for iSF = 1:size(theConeMosaicSTFresponses,2)
            
            r = theConeMosaicSTFresponses(7,iSF,:,:);
            
            c = squeeze(bsxfun(@times, bsxfun(@minus,r,n),1./n));
size(c)
    figure()
    c =c(:,idx);
    [min(c(:)) max(c(:))]
            STF = max(c,[],1);
            imagesc(c, [-1 1])
colormap(gray(1024))
            max(STF(:))
            size(STF)
            figure(99); clf;
            lConeIndices = find(theMidgetRGCMosaic.inputConeMosaic.coneTypes(idx) == cMosaic.LCONE_ID)
            mConeIndices = find(theMidgetRGCMosaic.inputConeMosaic.coneTypes(idx) == cMosaic.MCONE_ID)
            plot(theMidgetRGCMosaic.inputConeMosaic.coneRFpositionsDegs(idx(lConeIndices),1), STF(lConeIndices), 'rs-')
            hold on;
            plot(theMidgetRGCMosaic.inputConeMosaic.coneRFpositionsDegs(idx(mConeIndices),1), STF(mConeIndices), 'gs-')
            pause
            drawnow;
        end

        % Retrieve the default retinal RF params
        retinalRFparams = theMidgetRGCMosaic.defaultRetinalRFparams;

        % Retrieve the default visual RF properties
        visualRFparams = theMidgetRGCMosaic.defaultVisualRFparams;

        % Optionally, modify certain retinalRFparams and/or visualRFparams
        
        theMidgetRGCMosaic.optimizeConeWiring(...
            inputConeMosaicSTFresponsesFileName, ...
            retinalRFparams, ...
            visualRFparams);
    end


end
