function ReidShapleyMequenceRFmapSingleRun(theMRGCMosaicResponsesFullFileName, ...
	surroundConnectedParamsStruct, smoothFittedRFmaps, visualizeSmoothingAndFitting)

	load(theMRGCMosaicResponsesFullFileName, ...
		'theMRGCMosaic', 'stimParams', 'RFmappingParamsStruct', ...
 		'spatialSupportDegs', ...
	    'mSequenceIndicatorFunctions', ...
 		'theNoiseFreeSpatioTemporalMRCMosaicResponse', 'theMRGCMosaicResponseTemporalSupportSeconds');

	% Compute RF maps for all RGCs
    theRFmaps = cell(theMRGCMosaic.rgcsNum, 1);
    mSequenceStimsNum = size(theNoiseFreeSpatioTemporalMRCMosaicResponse,2);

    parfor iRGC = 1:theMRGCMosaic.rgcsNum
        fprintf('Computing msequence RF for mRGC %d of %d\n', iRGC, theMRGCMosaic.rgcsNum);
        theRFmap = zeros(size(mSequenceIndicatorFunctions,2), size(mSequenceIndicatorFunctions,3));
        allMsequenceStimResponses = theNoiseFreeSpatioTemporalMRCMosaicResponse(1, :, iRGC);

        for iStim = 1:mSequenceStimsNum
            theStimulusFrame = double(squeeze(mSequenceIndicatorFunctions(iStim,:,:)));
            theRFmap = theRFmap + theStimulusFrame * squeeze(allMsequenceStimResponses(iStim));
        end
        theRFmap = theRFmap / mSequenceStimsNum;

        % Expand to full size
        theRFmaps{iRGC} = visualStimulusGenerator.expandFrame(theRFmap, stimParams.rfPixelRetinalPixelsWithin);

    end % parfor iRGC

    % Analyze the maps by fitting a Gaussian ellipsoid
    GaussianEllipsoidSamplesNum = 256;
    onlyFitPositiveRFmap = true;
    fixCentroid = true;
    if (smoothFittedRFmaps)
        [RFmapRows, RFmapCols] = size(theRFmaps{1});
        smoothingKernelSizePixels = RFmapCols / (2.0*sqrt(numel(theRFmaps)));
    else
        smoothingKernelSizePixels = [];
    end

    
    theFittedGaussianEllipsoids = RGCMosaicAnalyzer.compute.GaussianEllipsoidFitsToSpatialRFmaps(...
        theRFmaps, spatialSupportDegs, GaussianEllipsoidSamplesNum, onlyFitPositiveRFmap, ...
        fixCentroid, smoothingKernelSizePixels, visualizeSmoothingAndFitting);

    % Update the ellipsoids with sigmas specified in microns
    for iRGC = 1:numel(theFittedGaussianEllipsoids)
        figure(22); clf;
        theRFmap = theRFmaps{iRGC};
        theFittedEllipsoid = theFittedGaussianEllipsoids{iRGC};
        if (isempty(theFittedEllipsoid))
            fprintf(2,'*** No fit for cell %d of %d\n', iRGC, numel(theFittedGaussianEllipsoids));
            continue
        end

        ax = subplot(1,2,1);
        imagesc(ax,spatialSupportDegs, spatialSupportDegs, theRFmap);
        set(ax, 'CLim', max(abs(theRFmap(:)))*[-1 1] );
        colormap(ax,gray)
        axis(ax, 'image');
        ax = subplot(1,2,2);
        imagesc(ax,theFittedEllipsoid.xSupportDegs, theFittedEllipsoid.ySupportDegs, theFittedEllipsoid.rfMap);
        set(ax, 'CLim', max(abs(theRFmap(:)))*[-1 1] );
        colormap(ax,gray)
        axis(ax, 'image');

        % Update the fitted Gaussian ellipsoid by adding the minor & major sigmas in microns
        updatedEllipsoid = theFittedGaussianEllipsoids{iRGC};
        if (isempty(updatedEllipsoid))
            continue;
        end

        % Compute the corresponding sigmas in microns
        updatedEllipsoid.minorSigmaMicrons = 1e3 * theMRGCMosaic.angularEccDegsToLinearEccMMs(updatedEllipsoid.minorSigmaDegs);
        updatedEllipsoid.majorSigmaMicrons = 1e3 * theMRGCMosaic.angularEccDegsToLinearEccMMs(updatedEllipsoid.majorSigmaDegs);
        theFittedGaussianEllipsoids{iRGC} = updatedEllipsoid;
    end

    % Compute the temporal equivalent eccentricities in degs
    theTemporalEquivalentEccentricitiesDegs = theMRGCMosaic.temporalEquivalentEccentricityForEccXYDegs(theMRGCMosaic.rgcRFpositionsDegs);

    % Compute the corresponding eccentricities in MMs
    theTemporalEquivalentEccentricitiesMMs = theMRGCMosaic.angularEccDegsToLinearEccMMs(theTemporalEquivalentEccentricitiesDegs);

    % Generate theRFmapsFileName
    matFileName = strrep(theMRGCMosaicResponsesFullFileName, 'mSequenceResponses', 'demos/ReidShapleyAnalyses');
    matFileName = strrep(matFileName, '.mat', 'RF.mat');
    matFileName = strrep(matFileName, '_mRGCMosaic', '');
    
    fprintf('\nSaving computed mRGCRF RF maps and their Gaussian ellipsoid fits to %s ...', matFileName);
    save(matFileName, ...
            'RFmappingParamsStruct', ...
            'theTemporalEquivalentEccentricitiesDegs', ...
            'theTemporalEquivalentEccentricitiesMMs', ...
            'theRFmaps', ...
            'spatialSupportDegs', ...
            'theFittedGaussianEllipsoids', ...
            '-v7.3');
end
