function [RsToRcVarianceCK, intStoCsensVarianceCK, RsToRcVariance, intStoCsensVariance, ...
	radialTemporalEquivalentEccentricityDegsMosaic, RcDegsMosaic, RsToRcMosaic, intStoCsensMosaic, KcMosaic] = CronerAndKaplanSTFanalysis(...
	theMRGCMosaicSTFResponsesFullFileName, ...
	theCronerKaplanAnalysisFileName, ...
	aggregatePreviouslyAnalyzedRunsFromMultipleTargetVisualSTFs, ...
    aggregatePreviouslyAnalyzedRunsFromMultipleEccentricities, ...
	targetedSurroundPurityRange , ...
    targetedRadialEccentricityRange, ...
    targetedCenterConeNumerosityRange, ...
    targetedCenterPurityRange , ...
    recomputeAnalysis, ...
    pdfExportSubDir, ...
    varargin)

	% Parse optional input
	p = inputParser;
	p.addParameter('limitVisualizedCKdataToTheEccentricititesOfSyntheticCells', false, @islogical);
	p.addParameter('fixedOptimalOrientation', [], @(x)(ischar(x)||isnumeric(x)||(isnan(x))));
	p.addParameter('deltaThresholdForLimitingFittedSTFtoPrimaryPeak', [], @isnumeric);
	p.addParameter('onlyDisplayCronerKaplanData', ~false, @islogical);
	p.addParameter('onlyVisualizeFittedDOGmodels', false, @islogical);
    p.addParameter('maxVisualizedRsRcRatio', 20, @isscalar);
    p.addParameter('maxVisualizedIntStoCsensRatio', 1.4, @isscalar);
    p.addParameter('visualizeMosaic', false, @islogical);
	p.addParameter('visualizeFullAndMaximalExcursionSTF', false, @islogical);
	p.addParameter('visualizeSTFfits', false, @islogical);
	p.addParameter('visualizeModelFitting', false, @islogical);
	p.addParameter('visualizeSTFwithConeWeightsMap', false, @islogical);
	p.addParameter('retinalSpaceReferredSTFs', false, @islogical);
	p.addParameter('mRGCNonLinearityParams', [], @(x)(isempty(x))||(isstruct(x)));
	p.addParameter('customTemporalFrequencyAndContrast', [], @(x)(isempty(x))||(isstruct(x)));
	p.addParameter('showComponentLineWeightingFunctions', false, @islogical);
    p.addParameter('generatePRCtilesPlots', false, @islogical);

    p.parse(varargin{:});
	limitVisualizedCKdataToTheEccentricititesOfSyntheticCells = p.Results.limitVisualizedCKdataToTheEccentricititesOfSyntheticCells;
	deltaThresholdForLimitingFittedSTFtoPrimaryPeak = p.Results.deltaThresholdForLimitingFittedSTFtoPrimaryPeak;
	fixedOptimalOrientation = p.Results.fixedOptimalOrientation;

	onlyDisplayCronerKaplanData = p.Results.onlyDisplayCronerKaplanData;
	onlyVisualizeFittedDOGmodels = p.Results.onlyVisualizeFittedDOGmodels;

    maxVisualizedRsRcRatio = p.Results.maxVisualizedRsRcRatio;
    maxVisualizedIntStoCsensRatio = p.Results.maxVisualizedIntStoCsensRatio;

    visualizeMosaic = p.Results.visualizeMosaic;
	visualizeFullAndMaximalExcursionSTF = p.Results.visualizeFullAndMaximalExcursionSTF;
	visualizeSTFfits = p.Results.visualizeSTFfits;
	visualizeModelFitting = p.Results.visualizeModelFitting;
	visualizeSTFwithConeWeightsMap = p.Results.visualizeSTFwithConeWeightsMap;
	retinalSpaceReferredSTFs = p.Results.retinalSpaceReferredSTFs;
	mRGCNonLinearityParams = p.Results.mRGCNonLinearityParams;
	customTemporalFrequencyAndContrast = p.Results.customTemporalFrequencyAndContrast;
	showComponentLineWeightingFunctions = p.Results.showComponentLineWeightingFunctions;
    generatePRCtilesPlots = p.Results.generatePRCtilesPlots;

	RsToRcVarianceCK = [];
	intStoCsensVarianceCK = [];
	RsToRcVariance = [];
	intStoCsensVariance = [];

	radialTemporalEquivalentEccentricityDegsMosaic = [];
	RcDegsMosaic = [];
	RsToRcMosaic = [];
	intStoCsensMosaic = [];
	KcMosaic = [];


    % Where PDF files are exported
    pdfExportRootDir = RGCMosaicConstructor.filepathFor.rawFigurePDFsDir();

    radialTemporalEquivalentEccentricityDegsAllRuns = [];

	if (recomputeAnalysis)
		if (iscell(theCronerKaplanAnalysisFileName))
			error('Expecting a single file when recomputing the CronerKaplan analysis\n');
		end
		fprintf('Recomputing Croner&Kaplan analysis from data in\n%s\n', theMRGCMosaicSTFResponsesFullFileName);

		% Load the responses
		load(theMRGCMosaicSTFResponsesFullFileName, 'theMRGCMosaic', 'stimParams', ...
			'theNoiseFreeSpatioTemporalMRGCMosaicResponses2DSTF', 'theMRGCMosaicResponseTemporalSupportSeconds');

        if (visualizeMosaic)
		    hFig = figure(1);
            set(hFig, 'Position', [10 10 1000 1000]);
            ax = subplot('Position', [0.05 0.05 0.94 0.94]);
            theMRGCMosaic.visualize(...
                'figureHandle', hFig, ...
                'axesHandle', ax, ...
                'centerSubregionContourSamples', 64, ...
                'spatialSupportSamples', 1024, ...
                'identifiedConeApertureThetaSamples', 64, ...
                'plottedRFoutlineFaceColor',  [0 0.8 0], ...
    		    'plottedRFoutlineFaceAlpha', 0.25, ...
    		    'identifyInputCones', true, ...
    		    'inputConesAlpha', 0.8, ...
    		    'identifiedConeAperture', 'lightCollectingAreaCharacteristicDiameter', ... %'geometricArea', ...
    		    'contourGenerationMethod', 'ellipseFitToPooledConePositions'); %'ellipseFitBasedOnLocalSpacing');  % 'contourOfPooledConeApertureImage'
    
            % Compute mosaic pdf filename from the theMRGCMosaicSTFResponsesFullFileName
            idx = strfind(theMRGCMosaicSTFResponsesFullFileName, 'mRGCMosaicSTFresponses');
            theMRGCMosaicSTFResponsesFullFileName = theMRGCMosaicSTFResponsesFullFileName(idx:end);
            postFix = strrep(strrep(theMRGCMosaicSTFResponsesFullFileName, '.mat', ''), 'mRGCMosaicSTFresponses_', '');
            pdfFileName = sprintf('%s.pdf', postFix);
        
            % Generate the path if we need to
            theVisualizationPDFfilename = fullfile(pdfExportSubDir, pdfFileName);
            % Generate the path if we need to
            thePDFFullFileName  = RGCMosaicConstructor.filepathFor.augmentedPathWithSubdirs(...
                    pdfExportRootDir, theVisualizationPDFfilename, ...
                    'generateMissingSubDirs', true);
    
            NicePlot.exportFigToPDF(thePDFFullFileName, hFig, 300, 'beVerbose');
        end


		% Preallocate memory
		optimalOrientationDegsMosaic = zeros(1, theMRGCMosaic.rgcsNum);
		KcMosaic = zeros(1, theMRGCMosaic.rgcsNum);
		intStoCsensMosaic = zeros(1, theMRGCMosaic.rgcsNum);
		RsToRcMosaic = zeros(1, theMRGCMosaic.rgcsNum);
		RcDegsMosaic = zeros(1, theMRGCMosaic.rgcsNum);
		radialTemporalEquivalentEccentricityDegsMosaic = zeros(1, theMRGCMosaic.rgcsNum);
		spatialFrequencyCPDFullRange = cell(1, theMRGCMosaic.rgcsNum);
	    theSTFtoFitFullRangeMosaic = cell(1, theMRGCMosaic.rgcsNum);
		spatialFrequencyCPD = cell(1, theMRGCMosaic.rgcsNum);
		theSTFtoFitMosaic = cell(1, theMRGCMosaic.rgcsNum);
		theFittedSTFsMosaic = cell(1, theMRGCMosaic.rgcsNum);
		
        
        if isempty(fixedOptimalOrientation)
            fprintf(2,'STFs will be analyzed at the orientation for which the STF at half max extends to the highest SF\n');
        elseif isnan(fixedOptimalOrientation)
            fprintf(2,'STFs will be analyzed at random orientations\n');
        elseif ischar(fixedOptimalOrientation)
	        fprintf(2,'STFs analyzed at the orientation that results in %s', fixedOptimalOrientation);
        else
            fprintf(2,'STFs will be analyzed at a fixed orientation: %2.1f degs\n', fixedOptimalOrientation);
        end
        pause(1.0);

        multiStartsNum = 128;

	    if (visualizeModelFitting) || (visualizeFullAndMaximalExcursionSTF) || ...
	    	   (visualizeSTFfits) || (visualizeSTFwithConeWeightsMap)

	    	% Visualizing stuff, so do serial computation
			for iRGC = 1:theMRGCMosaic.rgcsNum

				if (retinalSpaceReferredSTFs)
					fittedSTFsPDFfilename = sprintf('RetinalSpaceFittedSTFs_MosaicEccDegs_%2.2f_RGC_%d.pdf', theMRGCMosaic.eccentricityDegs(1),iRGC);
					pdfFileName = sprintf('RetinalSpaceSTFandConePoolingWeightsMap_MosaicEccDegs_%2.2f_RGC_%d.pdf', theMRGCMosaic.eccentricityDegs(1),iRGC);
				else
					fittedSTFsPDFfilename = sprintf('VisualSpaceFittedSTFs_MosaicEccDegs_%2.2f_RGC_%d.pdf', theMRGCMosaic.eccentricityDegs(1),iRGC);
					pdfFileName = sprintf('VisualSpaceSTFandConePoolingWeightsMap_MosaicEccDegs_%2.2f_RGC_%d.pdf', theMRGCMosaic.eccentricityDegs(1),iRGC);
				end


				% Update pdfFilenames to include the mRGCNonLinearity type and custom stim parameters
				[fittedSTFsPDFfilename, pdfFileName] = updatePSFfileNames(...
					fittedSTFsPDFfilename, pdfFileName, mRGCNonLinearityParams, customTemporalFrequencyAndContrast);


                theVisualizationPDFfilename = fullfile(pdfExportSubDir, fittedSTFsPDFfilename);
                % Generate the path if we need to
                thePDFFullFileName = RGCMosaicConstructor.filepathFor.augmentedPathWithSubdirs(...
                    pdfExportRootDir, theVisualizationPDFfilename, ...
                    'generateMissingSubDirs', true);

				dataOut = fitTheSTF(theMRGCMosaic, iRGC, stimParams, ...
					theMRGCMosaicResponseTemporalSupportSeconds, ...
					squeeze(theNoiseFreeSpatioTemporalMRGCMosaicResponses2DSTF(:,:,:,iRGC)), ...
					retinalSpaceReferredSTFs, ...
					fixedOptimalOrientation, ...
					deltaThresholdForLimitingFittedSTFtoPrimaryPeak, ...
					multiStartsNum, ...
					visualizeModelFitting, ...
					visualizeFullAndMaximalExcursionSTF, ...
					visualizeSTFfits, ...
					thePDFFullFileName);

				% Retrieve the data
				optimalOrientationDegsMosaic(iRGC) = dataOut.theOptimalOrientation;
				KcMosaic(iRGC) = dataOut.Kc;
				intStoCsensMosaic(iRGC) = dataOut.intStoCsens;
				RsToRcMosaic(iRGC) = dataOut.RsToRc; 
				RcDegsMosaic(iRGC) = dataOut.RcDegs; 
				radialTemporalEquivalentEccentricityDegsMosaic(iRGC) = dataOut.radialTemporalEquivalentEccentricityDegs;
				spatialFrequencyCPDFullRange{iRGC} = dataOut.spatialFrequencyCPDFullRange;
			    theSTFtoFitFullRangeMosaic{iRGC} = dataOut.theSTFtoFitFullRange;
				spatialFrequencyCPD{iRGC} = dataOut.spatialFrequencyCPD; 
				theSTFtoFitMosaic{iRGC} = dataOut.theSTFtoFit;
				theFittedSTFsMosaic{iRGC} = dataOut.theFittedSTFs;

				if (visualizeSTFwithConeWeightsMap)
                    theVisualizationPDFfilename = fullfile(pdfExportSubDir, pdfFileName);
                    % Generate the path if we need to
                    thePDFFullFileName = RGCMosaicConstructor.filepathFor.augmentedPathWithSubdirs(...
                        pdfExportRootDir, theVisualizationPDFfilename, ...
                        'generateMissingSubDirs', true);

					visualizeSingleCellConeWeightsMapsAndSTF(theMRGCMosaic, iRGC, dataOut, ...
						thePDFFullFileName, ...
						showComponentLineWeightingFunctions);
				end
			end % for iRGC = 1:theMRGCMosaic.rgcsNum
		else
			% Parallel execution if we are not visualizing stuff
			parfor iRGC = 1:theMRGCMosaic.rgcsNum
                if (retinalSpaceReferredSTFs)
					fittedSTFsPDFfilename = sprintf('RetinalSpaceFittedSTFs_MosaicEccDegs_%2.2f_RGC_%d.pdf', theMRGCMosaic.eccentricityDegs(1),iRGC);
					pdfFileName = sprintf('RetinalSpaceSTFandConePoolingWeightsMap_MosaicEccDegs_%2.2f_RGC_%d.pdf', theMRGCMosaic.eccentricityDegs(1),iRGC);
				else
					fittedSTFsPDFfilename = sprintf('VisualSpaceFittedSTFs_MosaicEccDegs_%2.2f_RGC_%d.pdf', theMRGCMosaic.eccentricityDegs(1),iRGC);
					pdfFileName = sprintf('VisualSpaceSTFandConePoolingWeightsMap_MosaicEccDegs_%2.2f_RGC_%d.pdf', theMRGCMosaic.eccentricityDegs(1),iRGC);
                end
               
               % Update pdfFilenames to include the mRGCNonLinearity type
				[fittedSTFsPDFfilename, pdfFileName] = updatePSFfileNames(...
					fittedSTFsPDFfilename, pdfFileName, mRGCNonLinearityParams, customTemporalFrequencyAndContrast);

                theVisualizationPDFfilename = fullfile(pdfExportSubDir, fittedSTFsPDFfilename);

                % Generate the path if we need to
                thePDFFullFileName = RGCMosaicConstructor.filepathFor.augmentedPathWithSubdirs(...
                    pdfExportRootDir, theVisualizationPDFfilename, ...
                    'generateMissingSubDirs', true);

				dataOut = fitTheSTF(theMRGCMosaic, iRGC, stimParams, ...
					theMRGCMosaicResponseTemporalSupportSeconds, ...
					squeeze(theNoiseFreeSpatioTemporalMRGCMosaicResponses2DSTF(:,:,:,iRGC)), ...
					retinalSpaceReferredSTFs, ...
					fixedOptimalOrientation, ...
					deltaThresholdForLimitingFittedSTFtoPrimaryPeak, ...
					multiStartsNum, ...
					visualizeModelFitting, ...
					visualizeFullAndMaximalExcursionSTF, ...
					visualizeSTFfits, ...
					thePDFFullFileName);

				% Retrieve the data
				optimalOrientationDegsMosaic(iRGC) = dataOut.theOptimalOrientation;
				KcMosaic(iRGC) = dataOut.Kc;
				intStoCsensMosaic(iRGC) = dataOut.intStoCsens;
				RsToRcMosaic(iRGC) = dataOut.RsToRc; 
				RcDegsMosaic(iRGC) = dataOut.RcDegs; 
				radialTemporalEquivalentEccentricityDegsMosaic(iRGC) = dataOut.radialTemporalEquivalentEccentricityDegs;
				spatialFrequencyCPDFullRange{iRGC} = dataOut.spatialFrequencyCPDFullRange;
			    theSTFtoFitFullRangeMosaic{iRGC} = dataOut.theSTFtoFitFullRange;
				spatialFrequencyCPD{iRGC} = dataOut.spatialFrequencyCPD; 
				theSTFtoFitMosaic{iRGC} = dataOut.theSTFtoFit;
				theFittedSTFsMosaic{iRGC} = dataOut.theFittedSTFs;
			end % parfor iRGC = 1:numel(targetRGCindices)
		end


		% Save the targetRGCindices that were analyzed
		save(theCronerKaplanAnalysisFileName, ...
            'fixedOptimalOrientation', ...
			'radialTemporalEquivalentEccentricityDegsMosaic', ...
			'RcDegsMosaic', ...
			'RsToRcMosaic', ....
			'intStoCsensMosaic', ...
			'KcMosaic', ...
			'optimalOrientationDegsMosaic', ...
			'spatialFrequencyCPDFullRange', ...
			'theSTFtoFitFullRangeMosaic', ...
			'spatialFrequencyCPD', ...
			'theSTFtoFitMosaic', ...
			'theFittedSTFsMosaic' ...
			);
		fprintf('Saved computed CronerKaplan analysis to %s\n', theCronerKaplanAnalysisFileName);

		visualizedRGCindices = 1:theMRGCMosaic.rgcsNum;
		intStoCsensMosaicEqualNumerosity = intStoCsensMosaic;
	end % if (recomputeAnalysis)



	if (~recomputeAnalysis)

		if (visualizeSTFwithConeWeightsMap)
			if (iscell(theCronerKaplanAnalysisFileName))
				runsNum = numel(theCronerKaplanAnalysisFileName);
			else
				runsNum = 1;
			end

			for iRun = 1:runsNum
				if (iscell(theCronerKaplanAnalysisFileName))
					theCurrentCronerKaplanAnalysisFileName= theCronerKaplanAnalysisFileName{iRun};
					theCurrentMRGCMosaicSTFResponsesFullFileName = theMRGCMosaicSTFResponsesFullFileName{iRun};
				else
					theCurrentCronerKaplanAnalysisFileName = theCronerKaplanAnalysisFileName;
					theCurrentMRGCMosaicSTFResponsesFullFileName = theMRGCMosaicSTFResponsesFullFileName;
				end
				fprintf('Loading analyzed C&K data from %s\n', theCurrentCronerKaplanAnalysisFileName);
				dataOut = load(theCurrentCronerKaplanAnalysisFileName);
				load(theCurrentMRGCMosaicSTFResponsesFullFileName, 'theMRGCMosaic');

				for iRGC = 1:theMRGCMosaic.rgcsNum
					if (retinalSpaceReferredSTFs)
						pdfFileName = sprintf('RetinalSpaceSTFandConePoolingWeightsMap_MosaicEccDegs_%2.2f_RGC_%d.pdf', ...
							theMRGCMosaic.eccentricityDegs(1),iRGC);
					else
						pdfFileName = sprintf('VisualSpaceSTFandConePoolingWeightsMap_MosaicEccDegs_%2.2f_RGC_%d.pdf', ...
							theMRGCMosaic.eccentricityDegs(1),iRGC);
					end

					% Update pdfFilenames to include the mRGCNonLinearity type
					[~, pdfFileName] = updatePSFfileNames(...
						'', pdfFileName, mRGCNonLinearityParams, customTemporalFrequencyAndContrast);
                    theVisualizationPDFfilename = fullfile(pdfExportSubDir, pdfFileName);
                
                    % Generate the path if we need to
                    thePDFFullFileName = RGCMosaicConstructor.filepathFor.augmentedPathWithSubdirs(...
                        pdfExportRootDir, theVisualizationPDFfilename, ...
                        'generateMissingSubDirs', true);

                    % Get analyzed STF data for this cell
                    dataOutForThisCell.spatialFrequencyCPDFullRange = dataOut.spatialFrequencyCPDFullRange{iRGC};
                    dataOutForThisCell.theSTFtoFitFullRange = dataOut.theSTFtoFitFullRangeMosaic{iRGC};
                    dataOutForThisCell.spatialFrequencyCPD = dataOut.spatialFrequencyCPD{iRGC};
                    dataOutForThisCell.theSTFtoFit = dataOut.theSTFtoFitMosaic{iRGC};
                    dataOutForThisCell.theFittedSTFs = dataOut.theFittedSTFsMosaic{iRGC};

                    if (isfield(dataOut,'optimalOrientationDegsMosaic'))
                    	dataOutForThisCell.optimalOrientationDegs = dataOut.optimalOrientationDegsMosaic(iRGC);
                    end
                    dataOutForThisCell.Kc = dataOut.KcMosaic(iRGC);
					dataOutForThisCell.intStoCsens = dataOut.intStoCsensMosaic(iRGC);
					dataOutForThisCell.RsToRc = dataOut.RsToRcMosaic(iRGC);
					dataOutForThisCell.RcDegs = dataOut.RcDegsMosaic(iRGC);

					visualizeSingleCellConeWeightsMapsAndSTF(theMRGCMosaic, iRGC, dataOutForThisCell, ...
						thePDFFullFileName, showComponentLineWeightingFunctions);
				end % iRGC
			end % for iRun
		end


        thePercentages = [5 10 25 50 75 90 95];

		if (aggregatePreviouslyAnalyzedRunsFromMultipleEccentricities) || (onlyVisualizeFittedDOGmodels)

            if (iscell(theCronerKaplanAnalysisFileName))
				runsNum = numel(theCronerKaplanAnalysisFileName);
			else
				runsNum = 1;
            end

            % Aggregated across all runs
			radialTemporalEquivalentEccentricityDegsMosaicAllRuns = [];
			RcDegsMosaicAllRuns = [];
			RsToRcMosaicAllRuns = [];
			intStoCsensMosaicAllRuns = [];
			intStoCsensMosaicIndividualRuns = {};
			KcMosaicAllRuns = [];

            % PRctiles for each run
            radialTemporalEquivalentEccentricityDegsAllRuns = cell(1,numel(runsNum));
            RcDegsPrcTilesAllRuns = zeros(runsNum, numel(thePercentages));
			RsToRcPrcTilesAllRuns = zeros(runsNum, numel(thePercentages));
			intStoCsensPrcTilesAllRuns = zeros(runsNum, numel(thePercentages));
			KcMosaicPrcTilesRuns = zeros(runsNum, numel(thePercentages));

			cellsNumIndividualRuns = [];
			desiredFixedOptimalOrientation = fixedOptimalOrientation;


			for iRun = 1:runsNum
				if (iscell(theCronerKaplanAnalysisFileName))
					theCurrentCronerKaplanAnalysisFileName= theCronerKaplanAnalysisFileName{iRun};
					theCurrentMRGCMosaicSTFResponsesFullFileName = theMRGCMosaicSTFResponsesFullFileName{iRun};
				else
					theCurrentCronerKaplanAnalysisFileName = theCronerKaplanAnalysisFileName;
					theCurrentMRGCMosaicSTFResponsesFullFileName = theMRGCMosaicSTFResponsesFullFileName;
				end
				fprintf('Loading analyzed C&K data from %s\n', theCurrentCronerKaplanAnalysisFileName);
				dataOut = load(theCurrentCronerKaplanAnalysisFileName);

				% Unpack the data for all cells in the mosaic
				fixedOptimalOrientation = dataOut.fixedOptimalOrientation;
				radialTemporalEquivalentEccentricityDegsMosaic = dataOut.radialTemporalEquivalentEccentricityDegsMosaic;
				RcDegsMosaic = dataOut.RcDegsMosaic;
				RsToRcMosaic = dataOut.RsToRcMosaic;
				intStoCsensMosaic = dataOut.intStoCsensMosaic;
				KcMosaic = dataOut.KcMosaic;

				spatialFrequencyCPDFullRange = dataOut.spatialFrequencyCPDFullRange;
			    theSTFtoFitFullRangeMosaic = dataOut.theSTFtoFitFullRangeMosaic;
				spatialFrequencyCPD = dataOut.spatialFrequencyCPD; 
				theSTFtoFitMosaic = dataOut.theSTFtoFitMosaic;
				theFittedSTFsMosaic = dataOut.theFittedSTFsMosaic;

	            if isempty(fixedOptimalOrientation)
	                	if (isempty(desiredFixedOptimalOrientation))
			                fprintf(2,'STFs analyzed at the orientation for which the STF at half max extends to the highest SF ');
			            else
			            	error('>>>>> orientation mismatch in computed STFs')
			            end
	            elseif isnan(fixedOptimalOrientation)
	                	 if (isempty(desiredFixedOptimalOrientation))
	                	 	error('>>>>> orientation mismatch in computed STFs')
	                	 elseif (isnan(desiredFixedOptimalOrientation))
	                	 	fprintf(2,'STFs analyzed at a random orientations');
	                	 else
	                	 	error('>>>>> orientation mismatch in computed STFs')
	                	 end
	            elseif ischar(fixedOptimalOrientation)
	                   fprintf(2,'STFs analyzed at the orientation that results in %s', fixedOptimalOrientation);
	            else
	                	if (isempty(desiredFixedOptimalOrientation))
	                		error('orientation mismatch in computed STFs')
	                	elseif (desiredFixedOptimalOrientation)
	                		error('>>>>> orientation mismatch in computed STFs')
	                	else
	                		if (fixedOptimalOrientation == desiredFixedOptimalOrientation)
	                    		fprintf(2,'STFs analyzed at a fixed orientation: %2.1f degs', fixedOptimalOrientation);
	                    	else
	                    		error('>>>> orientation mismatch in computed STFs')
	                    	end
	                    end
	            end

				load(theCurrentMRGCMosaicSTFResponsesFullFileName, 'theMRGCMosaic');

				% Select rgcIndices to visualize
				allRGCindices = 1:theMRGCMosaic.rgcsNum;
				targetRGCindices = theMRGCMosaic.indicesOfRGCsWithinTargetedPropertyRanges( ...
					targetedCenterConeNumerosityRange, ...
					targetedSurroundPurityRange, ...
					targetedRadialEccentricityRange, ...
					targetedCenterPurityRange);

				visualizedRGCindices = [];
				for i = 1:numel(targetRGCindices)
					idx = find(allRGCindices == targetRGCindices(i));
				    if (~isempty(idx))
						visualizedRGCindices(i) = allRGCindices(idx);
					end
				end

				if (onlyVisualizeFittedDOGmodels) && (~isempty(spatialFrequencyCPD))
					visualizeAllSTFfits(spatialFrequencyCPDFullRange, theSTFtoFitFullRangeMosaic, spatialFrequencyCPD, theSTFtoFitMosaic, theFittedSTFsMosaic, theMRGCMosaic, visualizedRGCindices);
				else
					fprintf('There are %d RGCs (out of a total of %d), that are within the targeted range\n', ...
						numel(visualizedRGCindices), theMRGCMosaic.rgcsNum);

					radialTemporalEquivalentEccentricityDegsMosaic = radialTemporalEquivalentEccentricityDegsMosaic(visualizedRGCindices);
					RcDegsMosaic = RcDegsMosaic(visualizedRGCindices);
					RsToRcMosaic = RsToRcMosaic(visualizedRGCindices);
					intStoCsensMosaic = intStoCsensMosaic(visualizedRGCindices);
					KcMosaic = KcMosaic(visualizedRGCindices);

                    % Compute the prctiles for this run
                    radialTemporalEquivalentEccentricityDegsAllRuns{iRun} = radialTemporalEquivalentEccentricityDegsMosaic(:);
                    RcDegsPrcTilesAllRuns(iRun,:) = prctile(RcDegsMosaic,thePercentages);
			        RsToRcPrcTilesAllRuns(iRun,:) = prctile(RsToRcMosaic,thePercentages);
			        intStoCsensPrcTilesAllRuns(iRun,:) = prctile(intStoCsensMosaic,thePercentages);
			        KcMosaicPrcTilesRuns(iRun,:) = prctile(KcMosaic,thePercentages);


					RsToRcVariance(iRun) = var(RsToRcMosaic);
					intStoCsensVariance(iRun) = var(intStoCsensMosaic);

					% Accumulate data over all runs
					radialTemporalEquivalentEccentricityDegsMosaicAllRuns = cat(1, ...
							radialTemporalEquivalentEccentricityDegsMosaicAllRuns(:), ...
							radialTemporalEquivalentEccentricityDegsMosaic(:));

					RcDegsMosaicAllRuns = cat(1, ...
							RcDegsMosaicAllRuns(:), RcDegsMosaic(:));

					RsToRcMosaicAllRuns = cat(1, ...
							RsToRcMosaicAllRuns(:), RsToRcMosaic(:));

					intStoCsensMosaicAllRuns = cat(1, ...
							intStoCsensMosaicAllRuns(:), intStoCsensMosaic(:));

					KcMosaicAllRuns  = cat(1, ...
							KcMosaicAllRuns(:), KcMosaic(:));

					intStoCsensMosaicIndividualRuns{iRun} = intStoCsensMosaic(:);
					cellsNumIndividualRuns(iRun) = numel(intStoCsensMosaic);
				end
			end % iRun

			if (~onlyVisualizeFittedDOGmodels)
				radialTemporalEquivalentEccentricityDegsMosaic = radialTemporalEquivalentEccentricityDegsMosaicAllRuns;
				RcDegsMosaic = RcDegsMosaicAllRuns;
				RsToRcMosaic = RsToRcMosaicAllRuns;
				KcMosaic = KcMosaicAllRuns;

				% intStoCsens: Case where we use all the data, despite the fact that different runs may have different # of cells
				intStoCsensMosaic = intStoCsensMosaicAllRuns;

				% intStoCsens: Case where use use equal # of cells across all eccentricities
				minCellsNum = min(cellsNumIndividualRuns);
				
				intStoCsensMosaicEqualNumerosity = [];
				for iRun = 1:numel(intStoCsensMosaicIndividualRuns)
					allData = intStoCsensMosaicIndividualRuns{iRun};
					randomizedCellIndices = randperm(numel(allData));
					equalNumerosityIntStoCsensMosaicData = allData(randomizedCellIndices(1:minCellsNum));
					intStoCsensMosaicEqualNumerosity = cat(1, ...
							intStoCsensMosaicEqualNumerosity(:), equalNumerosityIntStoCsensMosaicData(:));
				end

				fprintf('Equal numerosity intStoC has %d cells vs %d cells for the full dataset\n', numel(intStoCsensMosaicEqualNumerosity), numel(intStoCsensMosaic));
			end %  if (~onlyVisualizeFittedDOGmodels)
		else
			if (aggregatePreviouslyAnalyzedRunsFromMultipleTargetVisualSTFs)

				radialTemporalEquivalentEccentricityDegsMosaicAllRuns = [];
				RcDegsMosaicAllRuns = [];
				RsToRcMosaicAllRuns = [];
				intStoCsensMosaicAllRuns = [];
				KcMosaicAllRuns = [];
				for iRun = 1:numel(theCronerKaplanAnalysisFileName)
					fprintf('Loading analyzed C&K data from %s\n', theCronerKaplanAnalysisFileName{iRun});
					load(theCronerKaplanAnalysisFileName{iRun}, ...
                        		'fixedOptimalOrientation', ...
						'radialTemporalEquivalentEccentricityDegsMosaic', ...
						'RcDegsMosaic', ...
						'RsToRcMosaic', ....
						'intStoCsensMosaic', ...
						'KcMosaic');

	                    if isempty(fixedOptimalOrientation)
	                        fprintf(2,'STFs analyzed at the orientation for which the STF at half max extends to the highest SF ');
	                    elseif isnan(fixedOptimalOrientation)
	                         fprintf(2,'STFs analyzed at a random orientations');
	                    elseif ischar(fixedOptimalOrientation)
	                        fprintf(2,'STFs analyzed at the orientation that results in %s', fixedOptimalOrientation);
	                    else
	                        fprintf(2,'STFs analyzed at a fixed orientation: %2.1f degs', fixedOptimalOrientation);
	                    end

					radialTemporalEquivalentEccentricityDegsMosaicAllRuns = cat(1, ...
						radialTemporalEquivalentEccentricityDegsMosaicAllRuns(:), radialTemporalEquivalentEccentricityDegsMosaic(:));

					RcDegsMosaicAllRuns = cat(1, ...
						RcDegsMosaicAllRuns(:), RcDegsMosaic(:));
					RsToRcMosaicAllRuns = cat(1, ...
						RsToRcMosaicAllRuns(:), RsToRcMosaic(:));
					intStoCsensMosaicAllRuns = cat(1, ...
						intStoCsensMosaicAllRuns(:), intStoCsensMosaic(:));
					KcMosaicAllRuns  = cat(1, ...
						KcMosaicAllRuns(:), KcMosaic(:));

                     % Compute the prctiles for this run
                    radialTemporalEquivalentEccentricityDegsAllRuns{iRun} = radialTemporalEquivalentEccentricityDegsMosaic(:);
                    RcDegsPrcTilesAllRuns(iRun,:) = prctile(RcDegsMosaic,thePercentages);
			        RsToRcPrcTilesAllRuns(iRun,:) = prctile(RsToRcMosaic,thePercentages);
			        intStoCsensPrcTilesAllRuns(iRun,:) = prctile(intStoCsensMosaic,thePercentages);
			        KcMosaicPrcTilesRuns(iRun,:) = prctile(KcMosaic,thePercentages);


					RsToRcVariance(iRun) = var(RsToRcMosaic);
					intStoCsensVariance(iRun) = var(intStoCsensMosaic);

				end % iRun
				radialTemporalEquivalentEccentricityDegsMosaic = radialTemporalEquivalentEccentricityDegsMosaicAllRuns;
				RcDegsMosaic = RcDegsMosaicAllRuns;
				RsToRcMosaic = RsToRcMosaicAllRuns;
				intStoCsensMosaic = intStoCsensMosaicAllRuns;
				KcMosaic = KcMosaicAllRuns;
			else
				% No aggregated runs, just a single one
				load(theCronerKaplanAnalysisFileName, ...
                    'fixedOptimalOrientation' ,...
					'radialTemporalEquivalentEccentricityDegsMosaic', ...
					'RcDegsMosaic', ...
					'RsToRcMosaic', ....
					'intStoCsensMosaic', ...
					'KcMosaic');

	                if isempty(fixedOptimalOrientation)
	                        fprintf(2,'STFs analyzed at the orientation for which the STF at half max extends to the highest SF ');
	                elseif isnan(fixedOptimalOrientation)
	                        fprintf(2,'STFs analyzed at a random orientations');
	                elseif ischar(fixedOptimalOrientation)
	                        fprintf(2,'STFs analyzed at the orientation that results in %s', fixedOptimalOrientation);
	                else
	                        fprintf(2,'STFs analyzed at a fixed orientation: %2.1f degs', fixedOptimalOrientation);
	                end
				fprintf('Loaded previously computed CronerKaplan analysis from %s\n', theCronerKaplanAnalysisFileName);
			
				% Load the mosaic
				load(theMRGCMosaicSTFResponsesFullFileName, 'theMRGCMosaic');

				% Select rgcIndices to visualize
				allRGCindices = 1:theMRGCMosaic.rgcsNum;
				targetRGCindices = theMRGCMosaic.indicesOfRGCsWithinTargetedPropertyRanges( ...
					targetedCenterConeNumerosityRange, ...
					targetedSurroundPurityRange, ...
					targetedRadialEccentricityRange, ...
					targetedCenterPurityRange);

				visualizedRGCindices = [];
				for i = 1:numel(targetRGCindices)
					idx = find(allRGCindices == targetRGCindices(i));
				    if (~isempty(idx))
						visualizedRGCindices(i) = allRGCindices(idx);
					end
				end

				fprintf('There are %d RGCs (out of a total of %d), that are within the targeted range\n', ...
					numel(visualizedRGCindices), theMRGCMosaic.rgcsNum);

				radialTemporalEquivalentEccentricityDegsMosaic = radialTemporalEquivalentEccentricityDegsMosaic(visualizedRGCindices);
				RcDegsMosaic = RcDegsMosaic(visualizedRGCindices);
				RsToRcMosaic = RsToRcMosaic(visualizedRGCindices);
				intStoCsensMosaic = intStoCsensMosaic(visualizedRGCindices);
				KcMosaic = KcMosaic(visualizedRGCindices);

                % Compute the prctiles for this run
                radialTemporalEquivalentEccentricityDegsAllRuns{1} = radialTemporalEquivalentEccentricityDegsMosaic(:);
                RcDegsPrcTilesAllRuns(1,:) = prctile(RcDegsMosaic,thePercentages);
			    RsToRcPrcTilesAllRuns(1,:) = prctile(RsToRcMosaic,thePercentages);
			    intStoCsensPrcTilesAllRuns(1,:) = prctile(intStoCsensMosaic,thePercentages);
			    KcMosaicPrcTilesRuns(1,:) = prctile(KcMosaic,thePercentages);

				RsToRcVariance = var(RsToRcMosaic(~isoutlier(RsToRcMosaic, 'percentiles',[5 95])));
				intStoCsensVariance = var(intStoCsensMosaic(~isoutlier(intStoCsensMosaic, 'percentiles',[5 95])));
			end

			intStoCsensMosaicEqualNumerosity = intStoCsensMosaic;
		end % else
	end % if (~recomputeAnalysis)

	if (onlyVisualizeFittedDOGmodels)
		return;
	end

	% PLOT RESULTS
	if (onlyDisplayCronerKaplanData)
		radialTemporalEquivalentEccentricityDegsMosaic = [];
		RcDegsMosaic = [];
		RsToRcMosaic = [];
		intStoCsensMosaic = [];
	end

	% Rc (degs) as a function of temporal equivalent eccentricity
	figNo = 500; 

	[eccentricityDegsCKdata, RcDegsCKdata] = ...
		RGCmodels.CronerKaplan.digitizedData.parvoCenterRadiusAgainstEccentricity();

	minEccSupportDegs = 0.3;
	maxEccSupportDegs = 40;
    defaultEccSupportRange = [minEccSupportDegs maxEccSupportDegs];
    GodatEccSupportRange = [0.03 maxEccSupportDegs];
	maxRcDegs = 0.25;
	populationContourInsteadOfPointCloud = struct(...
		'xSupport', linspace(0, maxEccSupportDegs,160*8), ...
		'ySupport', linspace(0, maxRcDegs, 40), ...
		'percentileLevels', 5:10:95);


	% Data appearing in Fig. 10 (left panel) of Godat et al 2020
	% showing the Rc of the 13 mRGCs recorded by Godat as projected in visual space
	% using the M3 optics with a pupil size of 2.5 mm
	Godat2022data = [...
		6.384331573581E-3	9.568167101662E-3; ...
		6.793281606904E-3	9.569386943813E-3; ...
		7.228427041864E-3	9.176801456370E-3; ...
		3.542639075598E-2	7.307080153796E-3; ...
		5.678775577592E-2	6.724629356041E-3; ...
		5.678775577592E-2	6.118007802468E-3; ...
		7.650376142266E-2	8.567634560372E-3; ...
		7.189829714288E-2	7.551985289528E-3; ...
		7.370624358891E-2	6.588352247651E-3; ...
		7.940745894553E-2	6.589360196465E-3; ...
		8.554966620407E-2	6.801360489162E-3; ...
		8.242136620382E-2	5.688633576499E-3; ...
		9.216697631198E-2	5.342387669960E-3];
	eccentricityDegsGodat22data = Godat2022data(:,1);
	RcDegsGodat2022data = Godat2022data(:,2);


    % Comparison of Rc: synthetic vs macaque mRGCs 
    % (without the foveal data  of Godat et al 2022)
    pdfFileName = 'RcDegs.pdf'; 
    theVisualizationPDFfilename = fullfile(pdfExportSubDir, pdfFileName);
    % Generate the path if we need to
    thePDFFullFileName = RGCMosaicConstructor.filepathFor.augmentedPathWithSubdirs(...
                    pdfExportRootDir, theVisualizationPDFfilename, ...
                    'generateMissingSubDirs', true);

    exportVisualizationPDF = true;
	RGCMosaicAnalyzer.visualize.doubleScatterPlot(figNo, ...
		radialTemporalEquivalentEccentricityDegsMosaic, RcDegsMosaic, 'o', 6, [0.0 0.4 1], 0.02, 0.0, ...
		populationContourInsteadOfPointCloud, ...
		eccentricityDegsCKdata, RcDegsCKdata, [], ...
		's', 16, [0.8 0.8 0.8], 0.8, 1.0, [], ...
		defaultEccSupportRange, [0.01 0.03 0.1 0.3 1 3 10 30], [0 maxRcDegs], 0:0.05:0.25,...
		'log', 'linear', ...
		'temporal equivalent eccentricity (degs)', 'Rc (degs)', ...
		{sprintf('synthetic mRGCs, n=%d',numel(RcDegsMosaic)), 'macaque mRGCs (C&K''95)'}, ...
		exportVisualizationPDF, thePDFFullFileName);
 
    % Comparison of Rc: synthetic vs macaque mRGCs 
    % (with the foveal data  of Godat et al 2022)
    theVisualizationPDFfilename = fullfile(pdfExportSubDir, strrep(pdfFileName, '.pdf', 'withoutGodat2022.pdf'));
    thePDFFullFileName = RGCMosaicConstructor.filepathFor.augmentedPathWithSubdirs(...
                    pdfExportRootDir, theVisualizationPDFfilename, ...
                    'generateMissingSubDirs', false);

    figNo = figNo+1;
    RGCMosaicAnalyzer.visualize.doubleScatterPlot(figNo, ...
		radialTemporalEquivalentEccentricityDegsMosaic, RcDegsMosaic, 'o', 6, [0.0 0.4 1], 0.02, 0.0, ...
		populationContourInsteadOfPointCloud, ...
		eccentricityDegsCKdata, RcDegsCKdata, [], ...
		's', 16, [0.8 0.8 0.8], 0.8, 1.0, [], ...
		GodatEccSupportRange, [0.01 0.03 0.1 0.3 1 3 10 30], [0 maxRcDegs], 0:0.05:0.25,...
		'log', 'linear', ...
		'temporal equivalent eccentricity (degs)', 'Rc (degs)', ...
		{sprintf('synthetic mRGCs, n=%d',numel(RcDegsMosaic)), 'macaque mRGCs (C&K''95)'}, ...
		exportVisualizationPDF, thePDFFullFileName);

    theVisualizationPDFfilename = fullfile(pdfExportSubDir, strrep(pdfFileName, '.pdf', 'withGodat2022.pdf'));
    thePDFFullFileName = RGCMosaicConstructor.filepathFor.augmentedPathWithSubdirs(...
                    pdfExportRootDir, theVisualizationPDFfilename, ...
                    'generateMissingSubDirs', false);

    figNo = figNo+1;
	RGCMosaicAnalyzer.visualize.doubleScatterPlot(figNo, ...
		radialTemporalEquivalentEccentricityDegsMosaic, RcDegsMosaic, 'o', 6, [0.0 0.4 1], 0.02, 0.0, ...
		populationContourInsteadOfPointCloud, ...
		eccentricityDegsGodat22data, RcDegsGodat2022data, [], ...
		'o', 14, [1.0 0.2 0.7], 0.5, 0.5, [], ...
		GodatEccSupportRange, [0.01 0.03 0.1 0.3 1 3 10 30], [0 maxRcDegs], 0:0.05:0.25,...
		'log', 'linear', ...
		'temporal equivalent eccentricity (degs)', 'Rc (degs)', ...
		{sprintf('synthetic mRGCs, n=%d',numel(RcDegsMosaic)), 'macaque mRGCs (Godat et al. 2022)'}, ...
		exportVisualizationPDF, thePDFFullFileName);

    
    if (generatePRCtilesPlots)
        if (~isempty(radialTemporalEquivalentEccentricityDegsAllRuns))
            figNo = figNo+1;
            theVisualizationPDFfilename = fullfile(pdfExportSubDir, strrep(pdfFileName, '.pdf', 'PrcTiles.pdf'));
            thePDFFullFileName = RGCMosaicConstructor.filepathFor.augmentedPathWithSubdirs(...
                        pdfExportRootDir, theVisualizationPDFfilename, ...
                        'generateMissingSubDirs', false);
    
            RGCMosaicAnalyzer.visualize.scatterPrcTilesPlot(figNo, ...
		        radialTemporalEquivalentEccentricityDegsAllRuns, RcDegsPrcTilesAllRuns, 'o', 14, [0.0 0.4 1], 0.3, 1.0, ...
		        eccentricityDegsCKdata, RcDegsCKdata, [], ...
		        's', 16, [0.8 0.8 0.8], 0.8, 1.0, ...
		        defaultEccSupportRange, [0.01 0.03 0.1 0.3 1 3 10 30], [0 maxRcDegs], 0:0.05:0.25,...
		        'log', 'linear', ...
		        'temporal equivalent eccentricity (degs)', 'Rc (degs)', ...
		        {sprintf('synthetic, n=%d',numel(RcDegsMosaic)), 'Croner & Kaplan'}, ...
		        exportVisualizationPDF, thePDFFullFileName);
        end
    end


	% Rs/Rc ratio as a function of temporal equivalent eccentricity
	[eccentricityDegsCKdata, RcToRSCKdata] = ...
		RGCmodels.CronerKaplan.digitizedData.parvoCenterSurroundRadiusRatioAgainstEccentricity();
	RsRcCKdata = 1./RcToRSCKdata;

	populationContourInsteadOfPointCloud.ySupport = linspace(0, maxVisualizedRsRcRatio, 50);
	RsRcCKtarget = RGCmodels.CronerKaplan.constants.surroundToCenterRcRatio * ones(size(eccentricityDegsCKdata));

    pdfFileName = 'RsRcRatios.pdf';
    theVisualizationPDFfilename = fullfile(pdfExportSubDir, pdfFileName);
    thePDFFullFileName = RGCMosaicConstructor.filepathFor.augmentedPathWithSubdirs(...
                    pdfExportRootDir, theVisualizationPDFfilename, ...
                    'generateMissingSubDirs', false);

    figNo = figNo+1; 
	RGCMosaicAnalyzer.visualize.doubleScatterPlot(figNo, ...
		radialTemporalEquivalentEccentricityDegsMosaic, RsToRcMosaic, 'o', 6, [0.0 0.4 1], 0.02, 0.0, ...
		populationContourInsteadOfPointCloud, ...
		eccentricityDegsCKdata, RsRcCKdata, RsRcCKtarget, ...
		's', 16, [0.8 0.8 0.8], 0.8, 1.0, [], ...
		defaultEccSupportRange, [0.01 0.03 0.1 0.3 1 3 10 30], [0 maxVisualizedRsRcRatio], 0:2:maxVisualizedRsRcRatio, ...
		'log', 'linear', ...
		'temporal equivalent eccentricity (degs)', 'Rs/Rc', ...
		{sprintf('synthetic, n=%d',numel(RsToRcMosaic)), 'Croner & Kaplan'}, ...
		exportVisualizationPDF, pdfFileName);

    if (generatePRCtilesPlots)
        if (~isempty(radialTemporalEquivalentEccentricityDegsAllRuns))
            figNo = figNo+1;
    
            theVisualizationPDFfilename = fullfile(pdfExportSubDir, strrep(pdfFileName, '.pdf', 'PrcTiles.pdf'));
            thePDFFullFileName = RGCMosaicConstructor.filepathFor.augmentedPathWithSubdirs(...
                        pdfExportRootDir, theVisualizationPDFfilename, ...
                        'generateMissingSubDirs', false);
    
            RGCMosaicAnalyzer.visualize.scatterPrcTilesPlot(figNo, ...
		        radialTemporalEquivalentEccentricityDegsAllRuns, RsToRcPrcTilesAllRuns, 'o', 14, [0.0 0.4 1], 0.3, 1, ...
		        eccentricityDegsCKdata, RsRcCKdata, RsRcCKtarget, ...
		        's', 16, [0.8 0.8 0.8], 0.8, 1.0, ...
		        defaultEccSupportRange, [0.01 0.03 0.1 0.3 1 3 10 30], [0 maxVisualizedRsRcRatio], 0:2:maxVisualizedRsRcRatio,...
		        'log', 'linear', ...
		        'temporal equivalent eccentricity (degs)', 'Rs/Rc', ...
		        {sprintf('synthetic, n=%d',numel(RcDegsMosaic)), 'Croner & Kaplan'}, ...
		        exportVisualizationPDF, thePDFFullFileName);
        end
    end

	% Int S/C ratio as a function of temporal equivalent eccentricity
	[eccentricityDegsCKdata, intSCratioCKdata] = ...
            RGCmodels.CronerKaplan.digitizedData.parvoSurroundCenterIntSensisitivityRatioAgainstEccentricity();

    RsToRcVarianceCK = var(RsRcCKdata(~isoutlier(RsRcCKdata, 'percentiles',[5 95])));
	intStoCsensVarianceCK = var(intSCratioCKdata(~isoutlier(intSCratioCKdata, 'percentiles',[5 95])));


    if (limitVisualizedCKdataToTheEccentricititesOfSyntheticCells)
	    	minMosaicEcc = min(radialTemporalEquivalentEccentricityDegsMosaic);
	    	maxMosaicEcc = max(radialTemporalEquivalentEccentricityDegsMosaic);
	    	relevantCKdataIndices = find(eccentricityDegsCKdata>=minMosaicEcc & eccentricityDegsCKdata<=maxMosaicEcc);
	    	macaqueDataLegend =  sprintf('macaque (ecc:%2.1f-%2.1f)', minMosaicEcc, maxMosaicEcc);
	else
		relevantCKdataIndices = 1:numel(intSCratioCKdata);
		macaqueDataLegend = sprintf('Croner & Kaplan');
	end

	
	populationContourInsteadOfPointCloud.ySupport = linspace(0, maxVisualizedIntStoCsensRatio, 50);
	intStoCsensCKtarget = RGCmodels.CronerKaplan.constants.surroundToCenterIntegratedSensitivityRatioFromEccDegsForPcells(eccentricityDegsCKdata(relevantCKdataIndices));

    pdfFileName = 'intSCRatios.pdf';
    theVisualizationPDFfilename = fullfile(pdfExportSubDir, pdfFileName);
    thePDFFullFileName = RGCMosaicConstructor.filepathFor.augmentedPathWithSubdirs(...
                    pdfExportRootDir, theVisualizationPDFfilename, ...
                    'generateMissingSubDirs', false);

    figNo = figNo+1; 
    RGCMosaicAnalyzer.visualize.doubleScatterPlot(figNo, ...
		radialTemporalEquivalentEccentricityDegsMosaic, intStoCsensMosaic, 'o', 6, [0.0 0.4 1], 0.05, 0.0, ...
		populationContourInsteadOfPointCloud, ...
		eccentricityDegsCKdata(relevantCKdataIndices), intSCratioCKdata(relevantCKdataIndices), intStoCsensCKtarget, ...
		's', 16, [0.8 0.8 0.8], 0.8, 1.0, [], ...
		defaultEccSupportRange, [0.01 0.03 0.1 0.3 1 3 10 30], [0 maxVisualizedIntStoCsensRatio], 0:0.2:maxVisualizedIntStoCsensRatio, ...
		'log', 'linear', ...
		'temporal equivalent eccentricity (degs)', 'ISs/ISc', ...
		{sprintf('synthetic, n=%d',numel(intStoCsensMosaic)), macaqueDataLegend}, ...
		exportVisualizationPDF, thePDFFullFileName);

    if (~isempty(radialTemporalEquivalentEccentricityDegsAllRuns))
        figNo = figNo+1;

        theVisualizationPDFfilename = fullfile(pdfExportSubDir, strrep(pdfFileName, '.pdf', 'PrcTiles.pdf'));
        thePDFFullFileName = RGCMosaicConstructor.filepathFor.augmentedPathWithSubdirs(...
                    pdfExportRootDir, theVisualizationPDFfilename, ...
                    'generateMissingSubDirs', false);

        RGCMosaicAnalyzer.visualize.scatterPrcTilesPlot(figNo, ...
		    radialTemporalEquivalentEccentricityDegsAllRuns, intStoCsensPrcTilesAllRuns, 'o', 14, [0.0 0.4 1], 0.3, 1, ...
		    eccentricityDegsCKdata(relevantCKdataIndices), intSCratioCKdata(relevantCKdataIndices), intStoCsensCKtarget, ...
		    's', 16, [0.8 0.8 0.8], 0.8, 1.0, ...
		    defaultEccSupportRange, [0.01 0.03 0.1 0.3 1 3 10 30], [0 maxVisualizedIntStoCsensRatio], 0:0.2:maxVisualizedIntStoCsensRatio,...
		    'log', 'linear', ...
		    'temporal equivalent eccentricity (degs)', 'ISs/ISc', ...
		    {sprintf('synthetic, n=%d',numel(intStoCsensMosaic)), 'Croner & Kaplan'}, ...
		    exportVisualizationPDF, thePDFFullFileName);
    end

    % Histogram of Rs/Rc ratios
	

	theLegends = {sprintf('synthetic, n=%d',numel(RsToRcMosaic)), 'Croner & Kaplan'};
	theLegends = {};
    RsToRcBins = 0:2:maxVisualizedRsRcRatio;
    	
    pdfFileName = 'RsRcRatiosHistogram.pdf';
    theVisualizationPDFfilename = fullfile(pdfExportSubDir, pdfFileName);
    thePDFFullFileName = RGCMosaicConstructor.filepathFor.augmentedPathWithSubdirs(...
                    pdfExportRootDir, theVisualizationPDFfilename, ...
                    'generateMissingSubDirs', false);

    figNo = figNo+1; 
    RGCMosaicAnalyzer.visualize.doubleHistogramPlot(figNo, ...
	    	RsToRcMosaic, RsToRcBins, [0 0.4 1], 1.4, ...
	    	RsRcCKdata, RsToRcBins, [0.85 0.85 0.85], 2.0, ...
	    	[0 maxVisualizedRsRcRatio], 0:2.0:maxVisualizedRsRcRatio, [0 0.6], ...
	    	'Rs/Rc', 'frequency', ...
	    	theLegends, ...
	    	exportVisualizationPDF, thePDFFullFileName);


    % Histogram of intS/C ratios
	intStoCsensBins = 0:0.2:maxVisualizedIntStoCsensRatio;

	% For the intS/C ratios, which vary with eccentricity, we need to limit the # of cells
	% to be the same across all runs

	theLegends = {sprintf('synthetic, n=%d',numel(intStoCsensMosaic)), macaqueDataLegend};
	theLegends = {};

    pdfFileName = 'intSCratiosHistogram.pdf';
    theVisualizationPDFfilename = fullfile(pdfExportSubDir, pdfFileName);
    thePDFFullFileName = RGCMosaicConstructor.filepathFor.augmentedPathWithSubdirs(...
                    pdfExportRootDir, theVisualizationPDFfilename, ...
                    'generateMissingSubDirs', false);

    figNo = figNo+1; 
	RGCMosaicAnalyzer.visualize.doubleHistogramPlot(figNo, ...
		intStoCsensMosaicEqualNumerosity, intStoCsensBins, [0 0.4 1], 0.8*(intStoCsensBins(2)-intStoCsensBins(1)), ...
		intSCratioCKdata(relevantCKdataIndices), intStoCsensBins, [0.85 0.85 0.85], intStoCsensBins(2)-intStoCsensBins(1), ...
    	[intStoCsensBins(1) intStoCsensBins(end)], intStoCsensBins(1):intStoCsensBins(2)-intStoCsensBins(1):intStoCsensBins(end), ...
    	[0 0.6], ...
    	'ISs/ISc', 'frequency', ...
    	theLegends, ...
    	exportVisualizationPDF, thePDFFullFileName);
end


function dataOut = fitTheSTF(theMRGCMosaic, theRGCindex, stimParams, ...
	theMRGCMosaicResponseTemporalSupportSeconds, ...
	theNoiseFreeSpatioTemporalMRGCMosaicResponses2DSTF, ...
	retinalSpaceReferredSTFs, ...
	fixedOptimalOrientation, deltaThresholdForLimitingFittedSTFtoPrimaryPeak, multiStartsNum, ...
	visualizeModelFitting, visualizeFullAndMaximalExcursionSTF, visualizeSTFfits, thePDFfileName)

	fprintf('Analyzing RGC %d of %d...\n', theRGCindex, theMRGCMosaic.rgcsNum);

    if (visualizeModelFitting) || (visualizeSTFfits)
	    hFig = figure(100); clf;
	    set(hFig, 'Position', [10 10 1100 950], 'Color', [1 1 1]);
	    axFullSTF = subplot(2,2,1);
        axSTFslice = subplot(2,2,2);
        axSTFslicePortionFitted = subplot(2,2,3);
        axFittedSTFslice = subplot(2,2,4);
    else
        axFullSTF = [];
        axSTFslice = [];
        axSTFslicePortionFitted = [];
        axFittedSTFslice = [];
    end

	[theMaximalExcursionSTFamplitudeSpectrum, theOptimalOrientation, ...
	 theMaximalExcursionSTFphaseSpectrum, ...
	 theUnscaledMaximalExcursionSTFamplitudeSpectrum, ...
	 theFullSTFamplitudeSpectra, ...
	 theFullSTFphaseSpectra, hFigSinusoidalFits] = RGCMosaicConstructor.helper.simulateExperiment.maximalExcursionSTFfrom2DSTF(...
				stimParams.orientationDegs, stimParams.spatialFrequencyCPD, stimParams.spatialPhasesDegs, ...
				theMRGCMosaicResponseTemporalSupportSeconds, ...
				theNoiseFreeSpatioTemporalMRGCMosaicResponses2DSTF, ...
				'visualizeFullAndMaximalExcursionSTF', visualizeFullAndMaximalExcursionSTF, ...
				'visualizeSinusoidalFits', visualizeModelFitting, ...
				'fixedOptimalOrientation', fixedOptimalOrientation, ...
				'axFullSTF', axFullSTF, ...
    		    'axSTFslice',axSTFslice);



	if (visualizeModelFitting)
		hFig = figure(hFigSinusoidalFits);
		NicePlot.exportFigToPDF(strrep(thePDFfileName, '.pdf', 'SinusoidalFits.pdf'), hFigSinusoidalFits,  300, 'beVerbose');
	end


	optimaOrientationIndex = find(stimParams.orientationDegs == theOptimalOrientation);
	theSTFtoFitFullRange = theUnscaledMaximalExcursionSTFamplitudeSpectrum;
	spatialFrequencySupportCPDtoFitFullRange = stimParams.spatialFrequencyCPD;

	normFactor = max(theSTFtoFitFullRange(:));
	theSTFtoFitFullRange = theSTFtoFitFullRange / normFactor;
	
	temporalEquivalentXYeccDegs = ...
		theMRGCMosaic.temporalEquivalentEccentricityForEccXYDegs(theMRGCMosaic.rgcRFpositionsDegs(theRGCindex,:));
	radialTemporalEquivalentEccentricityDegs = sqrt(sum(temporalEquivalentXYeccDegs.^2,2));

	intStoCsens = struct(...
		'initial', RGCmodels.CronerKaplan.constants.surroundToCenterIntegratedSensitivityRatioFromEccDegsForPcells(radialTemporalEquivalentEccentricityDegs), ...
		'low', RGCmodels.CronerKaplan.constants.surroundToCenterIntegratedSensitivityRatioFromEccDegsForPcells(0)*0.2, ...
		'high', RGCmodels.CronerKaplan.constants.surroundToCenterIntegratedSensitivityRatioFromEccDegsForPcells(30)*5);

	RsToRc = struct(...
		'initial', RGCmodels.CronerKaplan.constants.surroundToCenterRcRatio, ...
		'low', RGCmodels.CronerKaplan.constants.surroundToCenterRcRatio/10, ...
		'high', RGCmodels.CronerKaplan.constants.surroundToCenterRcRatio*10);

	RsDegs = RGCmodels.CronerKaplan.constants.surroundCharacteristicRadiusFromFitToPandMcells(radialTemporalEquivalentEccentricityDegs);
	RcDegsInitial = RsDegs/RGCmodels.CronerKaplan.constants.surroundToCenterRcRatio;

	RcDegs = struct(...
		'initial', RcDegsInitial, ...
   	'low', RcDegsInitial/100, ...
   	'high', RcDegsInitial*100);

	initialKc = RGCmodels.CronerKaplan.constants.centerPeakSensitivityFromCharacteristicRadiusDegsForPcells(RcDegs.initial);
	Kc = struct(...
		'initial', initialKc, ... 
   	'low', initialKc/1000, ...
   	'high', initialKc*1000);

	%                          Kc           intStoCsens             RsToRc           RcDegs    
    	DoGparams.initialValues = [Kc.initial   intStoCsens.initial    RsToRc.initial    RcDegs.initial];
    	DoGparams.lowerBounds   = [Kc.low       intStoCsens.low        RsToRc.low        RcDegs.low];
    	DoGparams.upperBounds   = [Kc.high      intStoCsens.high       RsToRc.high       RcDegs.high];
    	DoGparams.names         = {'Kc',        'intStoCsens',         'RsToRc',         'RcDegs'};
    	DoGparams.scaling       = {'log',       'linear',           'linear',         'linear'};

    	if (visualizeModelFitting)
    		figure(55); clf;
    		axDoGFitToCompositeSTF = subplot(1,1,1);
    	else
    		axDoGFitToCompositeSTF = [];
    	end

	if (~isempty(deltaThresholdForLimitingFittedSTFtoPrimaryPeak))

		% Only fit the part of the STF that contains the main peak, avoiding secondary peaks
		% which are ususally the STF of individual cones in the RF center
		if (retinalSpaceReferredSTFs)
			% In retinal-space referred STFs, detect portion to fit using peak prominance
			minPeakProminance = 1e-2;
		else
			minPeakProminance = [];
		end
	    	[spatialFrequencySupportCPDtoFit, theSTFtoFit] = RGCMosaicConstructor.helper.simulateExperiment.stfPortionToAnalyze(...
	    				spatialFrequencySupportCPDtoFitFullRange, theSTFtoFitFullRange, ...
	    				deltaThresholdForLimitingFittedSTFtoPrimaryPeak, minPeakProminance);

	    	if (visualizeSTFfits)
	    		plot(axSTFslicePortionFitted, spatialFrequencySupportCPDtoFitFullRange, theSTFtoFitFullRange, 'ks-', 'LineWidth', 1.5)
        		set(axSTFslicePortionFitted, 'XScale', 'log', 'XLim', [0.1 100], 'YLim', [0 1]);
        		hold(axSTFslicePortionFitted, 'on')
        		plot(axSTFslicePortionFitted, spatialFrequencySupportCPDtoFit, theSTFtoFit, 'ro-', 'LineWidth', 1.5, 'MarkerSize', 15);
        		set(axSTFslicePortionFitted, 'FontSize', 16);
        		title(axSTFslicePortionFitted, 'STF slice portion to be fitted');
		end
	end

	[DoGparams, theFittedSTF] = RGCMosaicConstructor.helper.fit.genericDifferentOfGaussiansToCompositeSTF(...
			DoGparams, spatialFrequencySupportCPDtoFit, theSTFtoFit, axDoGFitToCompositeSTF, normFactor, multiStartsNum);

	% Form dataOut struct
	dataOut.spatialFrequencyCPDFullRange = spatialFrequencySupportCPDtoFitFullRange;
	dataOut.theOptimalOrientation = theOptimalOrientation;
	dataOut.theSTFtoFitFullRange = theSTFtoFitFullRange*normFactor;
	dataOut.spatialFrequencyCPD = spatialFrequencySupportCPDtoFit;
	dataOut.theSTFtoFit = theSTFtoFit*normFactor;
	dataOut.theFittedSTFs = theFittedSTF;

	idx = find(strcmp(DoGparams.names,'Kc'));
	dataOut.Kc = DoGparams.finalValues(idx);

	idx = find(strcmp(DoGparams.names, 'intStoCsens'));
	dataOut.intStoCsens = DoGparams.finalValues(idx);

	idx = find(strcmp(DoGparams.names, 'RsToRc'));
	dataOut.RsToRc = DoGparams.finalValues(idx);

	idx = find(strcmp(DoGparams.names, 'RcDegs'));
	dataOut.RcDegs = DoGparams.finalValues(idx);

	dataOut.radialTemporalEquivalentEccentricityDegs = radialTemporalEquivalentEccentricityDegs;


	if (visualizeSTFfits)
		theTitle = sprintf('RGC %d/%d at (%2.2f,%2.2f)', theRGCindex, theMRGCMosaic.rgcsNum, ...
			theMRGCMosaic.rgcRFpositionsDegs(theRGCindex,1), theMRGCMosaic.rgcRFpositionsDegs(theRGCindex,2));
		visualizeTheSTFfit(axFittedSTFslice, dataOut.spatialFrequencyCPDFullRange, dataOut.theSTFtoFitFullRange, ...
			dataOut.spatialFrequencyCPD, dataOut.theSTFtoFit, ...
			dataOut.theFittedSTFs, theTitle);

    	NicePlot.exportFigToPDF(thePDFfileName,hFig,  300, 'beVerbose');
	end
end


function visualizeSingleCellConeWeightsMapsAndSTF(theMRGCMosaic, iRGC, dataOut, thePDFfileName, showComponents)
	hFig = figure(1000); clf;
	set(hFig, 'Position', [10 10 1700 650], 'Color', [1 1 1]);
	axConeWeightsMap = axes('Position', [0.07 0.1 0.27 0.9]);
	axConeWeightsLineWeightingFunctions = axes('Position', [0.39 0.1 0.27 0.9]);
	axFittedSTFslice = axes('Position', [0.72 0.1 0.27 0.9]);

    % Compute default limits and ticks
	theRGCpositionDegs = theMRGCMosaic.rgcRFpositionsDegs(iRGC,:);
	[scaleBarDegs, scaleBarMicrons, ...
	 spatialSupportTickSeparationArcMin, ...
	 spatialSupportCenterDegs, ...
		 domainVisualizationLimits, ...
		 domainVisualizationTicks, ...
		 domainVisualizationLimitsSingleRF, ...
		 domainVisualizationTicksSingleRF] = ...
		 	RGCMosaicAnalyzer.visualize.generateLimits(theMRGCMosaic, theRGCpositionDegs);

    scaleBarDegs = [];

    % Visualize the mosaic of mRGC RF centers
    % identifying cones that are pooled by the RF center mechanism with
    % a weight >= mRGCMosaic.sensitivityAtPointOfOverlap;
    % This representation is like the representation used in visualizing 
    % mosaics of RGCs in typical in-vitro experiments (e.g. by the Chichilnisky lab)
    minCenterConeWeight = mRGCMosaic.sensitivityAtPointOfOverlap;

    % Include surround cones whose pooling weights are >= 0.001
    minSurroundConeWeight = 0.001;

    % Surround weight specified relative to the center weight, so common
    % scale with the center weights
    minSurroundConeWeightRelativity = 'center';


    theTitle = sprintf('RGC %d/%d @(%2.2f,%2.2f); Wc > %1.3f; Ws > %1.4f', iRGC, theMRGCMosaic.rgcsNum, ...
		theRGCpositionDegs(1), theRGCpositionDegs(2), ...
        minCenterConeWeight, minSurroundConeWeight);



	figureFormat = PublicationReadyPlotLib.figureComponents('1x1 standard figure'); 
	figNo = [];
    [centerLineWeightingFunctions, surroundLineWeightingFunctions] = ...
     	RGCMosaicAnalyzer.visualize.singleRGCconePoolingMap(figNo, ...
                theMRGCMosaic, iRGC, '', ...
                'domainVisualizationLimits', domainVisualizationLimitsSingleRF, ...
                'domainVisualizationTicks', domainVisualizationTicksSingleRF, ...
                'fixedSpatialSupportTickSeparationArcMin', spatialSupportTickSeparationArcMin, ...
                'fixedScaleBarDegs', scaleBarDegs, ...
                'doNotLabelScaleBar', true, ...
                'noGrid', true, ...
                'minCenterConeWeight', minCenterConeWeight, ...
                'minSurroundConeWeight', minSurroundConeWeight, ...
                'minSurroundConeWeightRelativity', minSurroundConeWeightRelativity, ...
                'plotTitle', theTitle, ...
                'figureHandle', hFig, ...
                'axesHandle', axConeWeightsMap, ...
                'figureFormat', figureFormat, ...
                'renderLineWeightingFunctionPlots', false);

     compositeInsteadOfComponentConePoolingLineWrightingFunction = false;
     whichMeridian = 'horizontal';
     RGCMosaicConstructor.visualize.centerAndSurroundConePoolingLineWeightingFunctions(...
			'', [], [], ...
			spatialSupportCenterDegs, spatialSupportTickSeparationArcMin, ...
			centerLineWeightingFunctions, surroundLineWeightingFunctions, whichMeridian, ...
			'axesToRenderIn', axConeWeightsLineWeightingFunctions, ...
			'domainVisualizationLimits', domainVisualizationLimitsSingleRF(1:2), ...
			'domainVisualizationTicks', domainVisualizationTicksSingleRF.x, ...
			'compositeInsteadOfComponent', compositeInsteadOfComponentConePoolingLineWrightingFunction);

    m1 = max([...
    	max(centerLineWeightingFunctions.xProfile.amplitude(:))
    	max(centerLineWeightingFunctions.yProfile.amplitude(:))
    	]);


    % Compute the surround line weighting function as derived from the DoG model's fit
	ks = dataOut.Kc * dataOut.intStoCsens/(dataOut.RsToRc^2);
	rsDegs = dataOut.RsToRc * dataOut.RcDegs;
	deMeanedSpatialSupportDegs = centerLineWeightingFunctions.xProfile.spatialSupportDegs-mean(centerLineWeightingFunctions.xProfile.spatialSupportDegs);
	[theCenterPointWeightingFunction, theSurroundPointWeightingFunction, theCenterLineWeightingFunction, theSurroundLineWeightingFunction, ...
	 theSTF, theSpatialSupport] = RGCmodels.EnrothCugellRobson.pointAndLineWeightingFunction(...
		deMeanedSpatialSupportDegs, ...
		[], dataOut.Kc, ks, dataOut.RcDegs, rsDegs, dataOut.spatialFrequencyCPDFullRange, false);


	% Normalize to same max
	m2 = max(theCenterLineWeightingFunction(:));
	theCenterLineWeightingFunction = theCenterLineWeightingFunction / m2 * m1;
	theSurroundLineWeightingFunction = theSurroundLineWeightingFunction / m2 * m1;


	hold(axConeWeightsLineWeightingFunctions, 'on')


	idx = find(theSurroundPointWeightingFunction >= 0.001*max(theSurroundPointWeightingFunction(:)));
	[rows,cols] = ind2sub([numel(deMeanedSpatialSupportDegs) numel(deMeanedSpatialSupportDegs)], idx);
	minX = min(deMeanedSpatialSupportDegs(rows));
	maxX = max(deMeanedSpatialSupportDegs(rows));
	idxSurround = find(...
		(deMeanedSpatialSupportDegs >= minX) & ...
		(deMeanedSpatialSupportDegs <= maxX));






	idx = find(theCenterPointWeightingFunction >= 0.0005*max(theCenterPointWeightingFunction(:)));
	[rows,cols] = ind2sub([numel(deMeanedSpatialSupportDegs) numel(deMeanedSpatialSupportDegs)], idx);

	minX = min(deMeanedSpatialSupportDegs(rows));
	maxX = max(deMeanedSpatialSupportDegs(rows));
	idx = find(...
		(deMeanedSpatialSupportDegs >= minX) & ...
		(deMeanedSpatialSupportDegs <= maxX));


	if (~showComponents)
		plot(axConeWeightsLineWeightingFunctions, ...
			surroundLineWeightingFunctions.xProfile.spatialSupportDegs(idxSurround), theCenterLineWeightingFunction(idxSurround)-theSurroundLineWeightingFunction(idxSurround), '-', 'Color', [0 0 0], 'LineWidth',4);

		plot(axConeWeightsLineWeightingFunctions, ...
			surroundLineWeightingFunctions.xProfile.spatialSupportDegs(idxSurround), theCenterLineWeightingFunction(idxSurround)-theSurroundLineWeightingFunction(idxSurround), '--', 'Color', [1 1 1], 'LineWidth', 2);

	else

		plot(axConeWeightsLineWeightingFunctions, ...
			centerLineWeightingFunctions.xProfile.spatialSupportDegs(idx), theCenterLineWeightingFunction(idx), '-', 'Color', [0 0 0], 'LineWidth',4);

		plot(axConeWeightsLineWeightingFunctions, ...
			centerLineWeightingFunctions.xProfile.spatialSupportDegs(idx), theCenterLineWeightingFunction(idx), '--', 'Color', [1 1 1], 'LineWidth', 2);

		plot(axConeWeightsLineWeightingFunctions, ...
			surroundLineWeightingFunctions.xProfile.spatialSupportDegs(idxSurround), -theSurroundLineWeightingFunction(idxSurround), '-', 'Color', [0 0 0], 'LineWidth',4);

		plot(axConeWeightsLineWeightingFunctions, ...
			surroundLineWeightingFunctions.xProfile.spatialSupportDegs(idxSurround), -theSurroundLineWeightingFunction(idxSurround), '--', 'Color', [1 1 1], 'LineWidth', 2);

	end



	hold(axConeWeightsLineWeightingFunctions, 'off')
	title(axConeWeightsLineWeightingFunctions, sprintf('DOG model Rs/Rc: %2.1f, Ks/Kc: %2.3fE-3', dataOut.RsToRc, ks/dataOut.Kc*1e3));

	theTitle = sprintf('fitted STF slice');
	visualizeTheSTFfit(axFittedSTFslice, dataOut.spatialFrequencyCPDFullRange, dataOut.theSTFtoFitFullRange, ...
		dataOut.spatialFrequencyCPD, dataOut.theSTFtoFit, ...
		dataOut.theFittedSTFs, theTitle);
	axis(axFittedSTFslice, 'square');

	PublicationReadyPlotLib.applyFormat(axFittedSTFslice,figureFormat);

	NicePlot.exportFigToPDF(thePDFfileName,hFig,  300, 'beVerbose');
	%RGCMosaicConstructor.helper.queryUserFor.unpausingExecution();
end 


function visualizeAllSTFfits(spatialFrequencyCPDFullRange, theSTFtoFitFullRangeMosaic, spatialFrequencyCPD, theSTFtoFitMosaic, theFittedSTFsMosaic, theMRGCMosaic, visualizedRGCindices)
	figure(2000); clf;
	ax = subplot(1,1,1);

	for iRGC = 1:numel(visualizedRGCindices)
		theRGCindex = visualizedRGCindices(iRGC);

		theSTFtoFitFullRange = theSTFtoFitFullRangeMosaic{theRGCindex};
		theSTFtoFit = theSTFtoFitMosaic{theRGCindex};
		theFittedSTF = theFittedSTFsMosaic{theRGCindex};

		theTitle = sprintf('RGC %d/%d at (%2.2f,%2.2f)', theRGCindex, numel(visualizedRGCindices), theMRGCMosaic.rgcRFpositionsDegs(theRGCindex,1), theMRGCMosaic.rgcRFpositionsDegs(theRGCindex,2));
		visualizeTheSTFfit(ax, spatialFrequencyCPDFullRange{theRGCindex}, theSTFtoFitFullRange, ...
			spatialFrequencyCPD{theRGCindex}, theSTFtoFit, theFittedSTF, theTitle);
	end % iRGC
end

function visualizeTheSTFfit(ax, spatialFrequencyCPDFullRange, theSTFtoFitFullRange, ...
		spatialFrequencyCPD, theSTFtoFit, theFittedSTF, theTitle)
	
	plot(ax, spatialFrequencyCPDFullRange, theSTFtoFitFullRange, 'ko', ...
		'MarkerFaceColor', [0.9 0.9 0.9], ...
		'MarkerEdgeColor', [0.75 0.75 0.75], ...
		'LineWidth', 1.5, ...
		'MarkerSize', 14);
	hold(ax, 'on')
	plot(ax, spatialFrequencyCPD, theSTFtoFit, 'ko', ...
		'MarkerFaceColor', [0.75 0.75 0.75], ...
		'MarkerEdgeColor', [0.5 0.5 0.5], ...
		'LineWidth', 2, ...
		'MarkerSize', 14);
	plot(ax, theFittedSTF.sfHiRes, theFittedSTF.compositeSTFHiRes, 'k-', 'LineWidth', 4.0);
	plot(ax, theFittedSTF.sfHiRes, theFittedSTF.centerSTFHiRes, 'k-', 'LineWidth', 3);
	plot(ax, theFittedSTF.sfHiRes, theFittedSTF.centerSTFHiRes, 'r-', 'Color', [1 0.5 0.5], 'LineWidth', 2);
	plot(ax, theFittedSTF.sfHiRes, theFittedSTF.surroundSTFHiRes, 'k-', 'LineWidth', 4.0);
	plot(ax, theFittedSTF.sfHiRes, theFittedSTF.surroundSTFHiRes, 'k-', 'Color', [166 191 210]/255,'LineWidth', 2);
	plot(ax, theFittedSTF.sfHiRes, theFittedSTF.compositeSTFHiRes, 'k-', 'LineWidth', 4.0);
	plot(ax, theFittedSTF.sfHiRes, theFittedSTF.compositeSTFHiRes, 'w--', 'LineWidth', 2.0);
	set(ax, 'XScale', 'log', 'XLim', [0.01 150], 'XTick', [0.03 0.1 0.3 1 3 10 30 100]);
	grid(ax, 'on');
	xlabel(ax, 'spatial frequency (c/deg)');
	title(ax, theTitle);
	drawnow
end


function [fittedSTFsPDFfilename, pdfFileName] = updatePSFfileNames(...
	fittedSTFsPDFfilename, pdfFileName, mRGCNonLinearityParams, customTemporalFrequencyAndContrast)


	if (~isempty(mRGCNonLinearityParams))
		% Update the fittedSTFsPDFfilename
		nonLinearityPreFix = 'FittedSTFs';
        switch (mRGCNonLinearityParams.type)
            case 'photocurrent'
                nonLinearityPostFix = 'FittedPhotocurrentSTFs';
            otherwise
                error('Unknown mRGCNonLinearityParams.type')
        end
        fittedSTFsPDFfilename = strrep(fittedSTFsPDFfilename, nonLinearityPreFix, nonLinearityPostFix);

        if (~isempty(customTemporalFrequencyAndContrast))
        	stimulusPreFix = 'FittedPhotocurrentSTFs';
        	stimulusPostFix = sprintf('%s%2.0f_x%1.0f_%02.0fHz', ...
                stimulusPreFix, ...
                customTemporalFrequencyAndContrast.totalContrast*100, ...
                customTemporalFrequencyAndContrast.backgroundLuminanceMultiplier, ...
                customTemporalFrequencyAndContrast.temporalFrequencyHz ...
                );
        	fittedSTFsPDFfilename = strrep(fittedSTFsPDFfilename, stimulusPreFix, stimulusPostFix);
        end

        % Update the pdfFileName
        nonLinearityPreFix = 'SpaceSTF';
        switch (mRGCNonLinearityParams.type)
            case 'photocurrent'
                nonLinearityPostFix = 'SpacePhotocurrentSTF';
            otherwise
                error('Unknown mRGCNonLinearityParams.type')
        end
        pdfFileName = strrep(pdfFileName, nonLinearityPreFix, nonLinearityPostFix);

        if (~isempty(customTemporalFrequencyAndContrast))
        	stimulusPreFix = 'SpacePhotocurrentSTF';
        	stimulusPostFix = sprintf('%s%2.0f_x%1.0f_%02.0fHz', ...
                stimulusPreFix, ...
                customTemporalFrequencyAndContrast.totalContrast*100, ...
                customTemporalFrequencyAndContrast.backgroundLuminanceMultiplier, ...
                customTemporalFrequencyAndContrast.temporalFrequencyHz ...
                );
        	pdfFileName = strrep(pdfFileName, stimulusPreFix, stimulusPostFix);
        end
	end

end
