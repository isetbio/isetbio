classdef MosaicPoolingOptimizer < handle

    properties (Constant)
        validChromaticityChoices = {'achromatic', 'Lcone isolating', 'Mcone isolating', 'Scone isolating'};
    end

    % Read-only properties 
    properties (SetAccess = private)
        % The RGC mosaic to be optimized
        theRGCMosaic;

        % The retinal RF model pooling params
        retinalRFmodelParams;

        % The nominal spatial sampling grid
        nominalSpatialSamplingGrid;

        % The sampling grids
        conesNumPooledByTheRFcenterGrid;
        targetRGCindicesWithLconeMajorityCenter;
        targetRGCindicesWithMconeMajorityCenter;
        visualSTFSurroundToCenterRcRatioGrid;
        visualSTFSurroundToCenterIntegratedSensitivityRatioGrid;
        
        % The computed L- and M-compute structs for all nodes in the sampling grid
        LconeComputeStructsList;
        MconeComputeStructsList;

        % The visual STF data of the input cone mosaic
        inputConeMosaicVisualSTFdata;

        % Number of multistarts for fitting the DoG model to the
        % model RGC STF
        multiStartsNumDoGFit;

        % Number of multistarts for fitting the surround model
        multiStartsNumRetinalPooling;

        % RMSE weights for the Rs/Rc and the S/CintSensitivity ratios
        rmseWeightForRsRcResidual;
        rmseWeightForSCintSensResidual;

        % Where to export PDFs
        exportedPDFFolder;
    end % Read only properties

    % Dependent properties
    properties (Dependent)
        % Number of grid nodes
        gridNodesNum;
    end


    properties (Constant)
        defaultRetinalRFmodelParams = struct(...
            'conePoolingModel', 'arbitraryCenterConeWeights_doubleExpH1cellIndex1SurroundWeights', ...  
            'H1parameterTolerance', 0.3, ...
            'centerConnectableConeTypes', [cMosaic.LCONE_ID cMosaic.MCONE_ID], ...   % cone types that can connect to the RF center
            'surroundConnectableConeTypes', [cMosaic.LCONE_ID cMosaic.MCONE_ID] ... % cone types that can connect to the RF surround
        );

        % H1 params for the 4 cells in Figure 6 of Packer & Dacey (2002)
        PackerDacey2002_H1params = struct(...
            'RnarrowToRwideRatios',  [152/515 170/718 115/902 221/1035], ...
            'NWvolumeRatios',        [1.0     0.8     0.3     0.2]);

        % Multiplier for surroundRc used to set the max radius for pooled surround
        % cones (modelConstants.maxSurroundSupportDegs)
        maxSurroundSupportFactor = 2.0;

        % Labels for available retinal cone pooling models
        ArbitraryCenter_DoubleExpH1cellIndex1Surround = 'arbitraryCenterConeWeights_doubleExpH1cellIndex1SurroundWeights';
        ArbitraryCenter_DoubleExpH1cellIndex2Surround = 'arbitraryCenterConeWeights_doubleExpH1cellIndex2SurroundWeights';
        ArbitraryCenter_DoubleExpH1cellIndex3Surround = 'arbitraryCenterConeWeights_doubleExpH1cellIndex3SurroundWeights';
        ArbitraryCenter_DoubleExpH1cellIndex4Surround = 'arbitraryCenterConeWeights_doubleExpH1cellIndex4SurroundWeights';

        % Attenuation factor at the high SF regime to determine optimal orientation
        deltaThresholdToDeclareLocalMinInSTF = 0.01;

    end % Constants

    % Public methods (class interface)
    methods
        % Constructor
        function obj = MosaicPoolingOptimizer(theRGCmosaic, ...
            varargin)

            % Parse input
            p = inputParser;
            p.addRequired('theRGCmosaic', @(x)(isa(x, 'mRGCMosaic')));
            p.addParameter('samplingScheme', 'hexagonal', @(x)(ismember(x, {'hexagonal', 'rectangular'})));
            p.addParameter('visualizeSamplingGrids', false, @islogical);
            p.addParameter('generateSamplingGrids', false, @islogical);
            p.addParameter('exportedPDFFolder', '', @ischar);
            p.parse(theRGCmosaic, varargin{:});

            obj.theRGCMosaic = p.Results.theRGCmosaic;
            obj.retinalRFmodelParams = obj.defaultRetinalRFmodelParams;
            obj.exportedPDFFolder = p.Results.exportedPDFFolder;

            if (p.Results.generateSamplingGrids)
                % Generate the nominal multifocal spatial sampling grid
                obj.generateNominalSpatialSamplingGrid(p.Results.samplingScheme);
                
                % Sampling positions must be at least 0.25 degs, or 10 x min RGC RF
                % spacing, whichever is larger
                minSpatialSamplingDegs = max([0.25 10*min(obj.theRGCMosaic.rgcRFspacingsDegs(:))]);

                % Generate the full multifocal sampling grids
                obj.generateSamplingGrids(minSpatialSamplingDegs);

                if (p.Results.visualizeSamplingGrids)
                    obj.visualizeSamplingGrids();
                end
            end

        end % Constructor

        % Method to visualize the sampling grids
        visualizeSamplingGrids(obj);

        % Method to compute the input cone mosaic visual STFs, either at a single
        % gridNode (using small stimulus patches centered on that node)
        % or across all nodes, using full field stimuli
        generateInputConeMosaicSTFresponses(obj, gridNodeIndex, ...
            stimSizeDegs, responsesFileName, varargin);

        % Method to compute the input cone mosaic visually projectedRcDegs by fitting a Gaussian to the
        % previously computed input cone mosaic visual STFs
        computeInputConeMosaicVisuallyProjectedRcDegs(obj, ...
            responsesFileName, varargin);

        % Method to load the computed input cone mosaic STF responses
        loadConeMosaicVisualSTFresponses(obj, responsesFileName);

        % Method to compute optimized RGC models (surround cone pooling weights)
        % for a grid node in the RGC mosaic
        compute(obj, gridNodeIndex, whichConeType, mosaicParams, opticsParams, ...
            coneMosaicSTFresponsesFileName, ...
            optimizedRGCpoolingObjectsFileName, ...
            varargin);

        % Method to inspect optimized RGC models for a grid node in the RGC mosaic
        inspect(obj, gridNodeIndex, opticsParams, ...
            optimizedRGCpoolingObjectsFileName, ...
            varargin);

        % Method to visualize the convergence of the optimized RGC model for a grid node
        animateModelConvergence(obj, gridNodeIndex, ...
            optimizedRGCpoolingObjectsFileName, ...
            varargin);

        % Method that incorporates the RGC pooling models optimized across all
        % nodes in the RGCmosaic and generates a compute-ready mRGCMosaic
        % with center and surround cone connections derived from these
        % optimized RGC pooling models
        generateComputeReadyMidgetRGCMosaic(obj, ...
            optimizedRGCpoolingObjectsFileName, ...
            computeReadyMosaicFilename, ...
            visualizeInterpolation);

       
        % Getter for dependent property gridNodesNum
        function val = get.gridNodesNum(obj)
            val = numel(obj.targetRGCindicesWithLconeMajorityCenter);
        end


    end % Public methods

    % Private methods
    methods (Access=private)
        % Method to generate the nominal spatial sampling grid
        generateNominalSpatialSamplingGrid(obj, samplingScheme);

        % Method to generate the full sampling grids
        generateSamplingGrids(obj, minSpatialSamplingDegs);

        % Method to return the center majority cone types
        [theCenterConeTypeWeights, theCenterConeTypeNum, theMajorityConeType, theCenterConeTypes] = ...
            centerConeTypeWeights(obj, theRGCindex);

      
        % Method to compute the RGCmodel STF and its DoG model fit params
        % based on its current cone pooling weights
        theSTFdata = rgcSTFfromPooledConeMosaicSTFresponses(obj, pooledConeIndicesAndWeights, visualRcDegs);

        % Method to opimize the surround cone pooling so as to achieve a
        % visual STF matching the targetSTF
        theRFcomputeStruct = optimizeSurroundConePooling(obj, theRGCindex, targetVisualSTFparams, ...
            mosaicParams, opticsParams, initialRetinalConePoolingParams, ...
            displayFittingProgress, exportedFittingProgressFolder, figNo, figTitle);

        % Optimization components
        [modelConstants, retinalConePoolingParams, visualRcDegs] = computeOptimizationComponents(...
            obj, theRGCindex, visualizeComponents, exportedFittingProgressFolder);
    end

    % Static methods
    methods (Static)

        % Method to generate the local dropbox path
        dropboxDir = localDropboxPath();

        % Fit a sinusoid to a response time series
        [theFittedResponse, fittedParams] = fitSinusoidToResponseTimeSeries(time, theResponse, stimulusTemporalFrequencyHz, timeHR);

        % Fit a Gaussian ellipsoid to a RFmap
        theFittedGaussianEllipsoid = fitGaussianEllipsoid(supportX, supportY, theRFmap, varargin)

        % Fit a Gaussian model to a subregion STF
        [Gparams, theFittedSTF] = fitGaussianToSubregionSTF(...
            spatialFrequencySupportCPD, theSubregionSTF, ...
            RcDegsInitialEstimate, rangeForRc, multiStartsNum)

        % Fit a DoG model to an STF
        [DoGparams, theFittedSTF] = fitDifferenceOfGaussiansToSTF(...
             spatialFrequencySupportCPD, theSTF, ...
             RcDegsInitialEstimate, rangeForRc, multiStartsNum);

        % Method to select STF portion to analyze (up to first local min)
        [sfSupport, theSTF] = stfPortionToAnalyze(sfSupport, theSTF)

        % Method to select the highest-extending STF (across a set of STFs
        % measured at different orientations)
        [theOptimalNormalizedSTF, theNormalizedSTFsAcrossAllOrientations, theHighestExtensionOrientation, ...
            theOptimalSTFMagnitudeSpectrum, theOptimalSTFphaseSpectrum] = optimalSTFfromResponsesToAllOrientationsAndSpatialFrequencies(...
                orientationsTested, spatialFrequenciesTested, ...
                theResponsesAcrossAllOrientationsAndSpatialFrequencies);


        % Method to compute the bandpass index of an STF as defined in:
        % Lee, Shapley, Hayken, and Sun. (2012) "Spatial distributions of cone inputs to
        % cells of the parvocellular pathway investigated with cone-isolating gratings", JOSA, 29, 2.
        BPI = STFbandpassIndex(theSpatialFrequencySupport, theSTF);

        % Method to convert cone pooling params to pooled cone indices and
        % weights for the double exponent H1 surround model 
        pooledConeIndicesAndWeights = conePoolingCoefficientsForArbitraryCenterDoubleExpSurround( ...
            modelConstants, conePoolingParamsVector);

        % Connectable cone indices to the surround
        [surroundConeIndices, surroundConeWeights, ...
          nonConnectableSurroundConeIndices, ...
          nonConnectableSurroundConeWeights] = connectableSurroundConeIndicesAndWeights(...
                surroundConeIndices, surroundConeWeights, modelConstants);

        % Method to compute triangulating nodes and weights for computing
        % pooling parameters from nearby RGC models
        [triangulatingModelRGCIndices, triangulatingModelRGCWeights, triangulatingGridNodeIndices] = ...
            triangulatingGridNodeIndicesAndWeights(theRGCposition, positionsOfModelRGCs, ...
            targetModelRGCindices, targetGridNodeIndices);

        % Method to visualize the optimization progress
        [hFig, figureFormat] = visualizeOptimizationProgress(figNo, figTitle, ...
            targetVisualSTFparams, opticsParams, theCurrentSTFdata, ...
            retinalConePoolingParams,  retinalConePoolingModel, ...
            pooledConeIndicesAndWeights, ...
            rmseSequence)

        % Method to visualize a model's parameter values & ranges
        visualizeFittedModelParametersAndRanges(ax, modelParams, modelName);

        % Method to visualize how we synthesize RF center/surround via
        % triangulating interpolation
        visualizeMultiSpectralMultiFocalRFgeneration(figNo, ...
            inputConeMosaic, theCurrentRGCindex, theCurrentRGCposition, ...
            coneType, ...
            coneIndicesAndWeightsForTriangulatingModels, ...
            triangulatingModelRGCpositions, ...
            triangulatingModelGridIndices, triangulatingModelWeights, ...
            centerConeIndicesForCurrentRGC, centerConeWeightsForCurrentRGC, ...
            surroundConeIndicesForCurrentRGC, surroundConeWeightsForCurrentRGC);

        % Method to visualize the retinal RF cone pooling map and visual
        % STF for a target RGC
        visualizeConePoolingRFmapAndVisualSTFforTargetRGC(...
            computeReadyMosaicFilename, mRGCMosaicSTFresponsesFilename, pdfFileName, ...
            targetRGCposition, targetCenterConesNum, targetCenterConeMajorityType, varargin);

        % Method to visualize the visual RF map computed via some way, such
        % as subspace mapping. The VisualRFmapStruct must contain the
        % following fields:
        % 'spatialSupportDegsX'   - vector
        % 'spatialSupportDegsY'   - vector
        % 'theRFmap'              - matrix 
        visualizeVisualRFmap(theVisualRFmapStruct, retinalRGCRFposDegs, theAxes, varargin);

        % Method to retrieve a subspace RF map for a target RGC
        retrieveAndVisualizeSubspaceRFmapForTargetRGC(...
            theComputeReadyMRGCmosaic, ...
            optimallyMappedSubspaceRFmapsFileName, ...
            targetRGCposition, targetCenterConesNum, targetCenterConeMajorityType, ...
            pdfFileName, varargin);

        % Method to visualize the visual RF maps for multiple target RGCs
        visualizeVisualRFmapsForMultipleTargetRGCs(...
            theComputeReadyMRGCmosaic, ...
            optimallyMappedSubspaceRFmapsFileName, ...
            pdfFileName, varargin);

     
        % Method to setup the parameters and the display for conducting an
        % STF mapping experiment.
        [stimParams, thePresentationDisplay] = setupSTFmappingExperiment(inputConeMosaic, ...
            sceneFOVdegs, retinalImageResolutionDegs, stimulusChromaticity);


        % Method to return cone contrasts and overall contrast for each
        % stimulus chromaticity
        [coneContrasts, contrast] = contrastForChromaticity(stimulusChromaticity);

        % Method to compute visualSTF responses of a cone mosaic under
        % some optics
        [theConeMosaicSTFresponses, theConeMosaicNullResponses] = ...
            computeConeMosaicSTFresponses(theConeMosaic, theOptics, ...
                                       thePresentationDisplay, ...
                                       stimParams, stimPositionDegs, ...
                                       useParfor, visualizeResponses);
        
        % Method to compute the visualSTF of the compute-ready
        % MRGCmosaic
        computeVisualSTFsOfComputeReadyMidgetRGCMosaic(...
            computeReadyMosaicFilename, ...
            coneMosaicSTFresponsesFilename, ...
            mRGCMosaicSTFresponsesFilename);

        % Method to fit the visualSTF of the compute-ready MRGCmosaic
        fitVisualSTFsOfComputeReadyMidgetRGCMosaic(...
                computeReadyMosaicFilename, ...
                mRGCMosaicSTFresponsesFilename);

        
        % Method to inspect the DoG model fits to measured STFs
        inspectDoGmodelFitsToMeasuredSTFs(computeReadyMosaicFilename, mRGCMosaicSTFresponsesFilename);

        % Method to visualize the fitted the visualSTF of the compute-ready MRGCmosaic
        visualizeFittedSTFsOfComputeReadyMidgetRGCMosaic(...
                computeReadyMosaicFilename, ...
                mRGCMosaicSTFresponsesFilename, ...
                pdfsDirectory, ...
                showZscoresInsteadOfData, visualizeFittedSTFsOfComputeReadyMidgetRGCMosaic, ...
                employTemporalEquivalentEccentricity);

        % Method to visualize cone mosaic STF responses
        visualizeConeMosaicSTFresponses(mRGCMosaicFileName, coneMosaicSTFresponsesFileName, varargin);

        % Method to visualize the vLambda weighted PSF for a compute-ready
        % mRGCmosaic given its opticsParams
        visualizeVlambdaWeightedPSF(theComputeReadyMRGCmosaic, opticsParams, varargin);

        % Method to contrast STFs across different optics/stimulus
        % chromaticities
        contrastVisualSTFsAcrossDifferentChromaticities(...
            exampleLconeCenterRGCposition, exampleMconeCenterRGCposition, ...
            theComputeReadyMRGCmosaic, ...
            mRGCMosaicAchromaticSTFresponsesFileName, ...
            mRGCMosaicLconeIsolatingSTFresponsesFileName, ...
            mRGCMosaicMconeIsolatingSTFresponsesFileName, ...
            opticsParams, rawFiguresRoot, scaledFiguresRoot, ...
            exportScaledFigureVersionForManuscript, ...
            varargin);


        % Method to contrast m-sequence RFs across different optics/stimulus
        % chromaticities
        contrastMSequenceRFsAcrossDifferentChromaticities(...
            exampleLconeCenterRGCposition, exampleMconeCenterRGCposition, ...
            theComputeReadyMRGCmosaic, ...
            theOptimallyMappedAchromaticRFmapsFileName, ...
            theOptimallyMappedLconeIsolatingRFmapsFileName, ...
            theOptimallyMappedMconeIsolatingRFmapsFileName, ...
            coneFundamentalsOptimizedForStimPosition, ...
            rfPixelsAcross, opticsParams, rawFiguresRoot, scaledFiguresRoot, ...
            exportScaledFigureVersionForManuscript, ...
            varargin);

        % Method to compute the surround cone mix for mRGCs
        [theCenterMajorityConeType, netCenterLconeWeight, netCenterMconeWeight, ...
         netSurroundLconeWeight, netSurroundMconeWeight, surroundConeMix, centerConeMix] = analyzeCenterSurroundConeMix(...
            theMRGCmosaic, theRGCindex, performSurroundAnalysisForConesExclusiveToTheSurround);

         % Method to determine the indices of mRGCs whose surround cone mix is within the target range
        mRGCindicesToVisualizeSTFsAcrossChromaticities = findMRGCindicesWithDesiredSurroundConeMix(...
            theComputeReadyMRGCmosaic, ...
            targetRangeForSurroundConeMix, ...
            targetRangeForCenterConeMix, ...
            performSurroundAnalysisForConesExclusiveToTheSurround, ...
            rawFiguresRoot, scaledFiguresRoot, exportScaledFigureVersionForManuscript);


        % Method to generate Vlambda weigted PSF for the mRGCmosaic
        thePSFData = generateVlambdaWeightedPSFData(theMidgetRGCMosaic, opticsParams)

        % Method to compute an optimal retinal image resolution for
        % stimulus generation (pixel size) based on which RGCs are to be
        % mapped (either by STF or Subspsace mapping)
        retinalImageResolutionDegs = retinalResolutionFromConeApertureDiameter(theMidgetRGCMosaic, targetRGCindices);

        % Method to derive cone fundamentals for a target region within the
        % input cone mosaic
        customConeFundamentals = coneFundamentalsAtTargetPositionWithinConeMosaic(...
            theConeMosaic, theOptics, targetRegionPositionDegs, targetRegionSizeDegs, maxConesNumForAveraging);

        % Method to compute cone PSFs at the stimulus position
        [theLconePSF, theMconePSF, theSconePSF, spatialSupportDegs] = ...
            computeConePSFs(theRGCMosaic, opticsToEmploy, stimSizeDegs, stimXYpositionDegs, varargin);

        % Method to compute the indices of optimally mapped RGCs within an
        % MRGC mosaic for a RF map conducted at stimPositionDegs
        indicesOfOptimallyMappedRGCs = indicesOfOptimallyMappedRGCsAtThisPosition(theMRGCmosaic, ...
            stimPositionDegs, stimSizeDegs);

        % Method to setup the parameters and the display for conducting an
        % Subspace RF mapping experiment.
        [stimParams, thePresentationDisplay] = setupSubspaceRFmappingExperiment(...
            wavelengthSupport, stimSizeDegs, ...
            optimalRetinalPixelSizeDegs, maxSFLimit, stimulusChromaticity);

        % Method to setup the parameters and the display for conducting an
        % M-sequence RF mapping experiment.
        [stimParams, thePresentationDisplay] = setupMSequenceRFmappingExperiment(...
              wavelengthSupport, stimSizeDegs, rfPixelsAcross, ...
              optimalRetinalPixelSizeDegs, ...
              ternaryInsteadOfBinaryMsequence, mSequenceBitLength, ...
              stimulusChromaticity);


        % Method to compute visual RFs of the computeReadyMosaic (using the
        % subspace RF mapping method)
        computeVisualRFsUsingSubSpaceMapping(...
            theComputeReadyMRGCmosaic, opticsToEmploy, ...
            stimSizeDegs, stimPositionDegs, ...
            maxSFLimit, maxSFToBeAnalyzed, rfMappingPixelMagnificationFactor, ...
            stimChromaticity, coneFundamentalsOptimizedForStimPosition, ...
            coneMosaicResponsesFileName, ...
            mRGCMosaicResponsesFileName, ...
            optimallyMappedRFmapsFileName, ...
            reComputeInputConeMosaicResponses, ...
            reComputeMRGCMosaicResponses, ...
            reComputeRFs, ...
            visualizeOptimallyMappedRFmapLocations, ...
            varargin);

        % Method to compute visual RFs of the computeReadyMosaic (using the
        % m-sequence RF mapping method)
        computeVisualRFsUsingMSequenceMapping(...
            theComputeReadyMRGCmosaic, opticsToEmploy, ...
            stimSizeDegs, stimPositionDegs, rfPixelsAcross, ...
            ternaryInsteadOfBinaryMsequence, mSequenceBitLength, ...
            stimChromaticity, coneFundamentalsOptimizedForStimPosition, ...
            coneMosaicResponsesFileName, ...
            mRGCMosaicResponsesFileName, ...
            optimallyMappedRFmapsFileName, ...
            reComputeInputConeMosaicResponses, ...
            reComputeMRGCMosaicResponses, ...
            reComputeRFs, ...
            visualizeOptimallyMappedRFmapLocations, ...
            varargin);

        % Method to visualize the visual RF centers of all RGCs in a
        % compute ready MRGC mosaic
        visualizeVisualRFcentersOfComputeReadyMidgetRGCMosaic(...
            theComputeReadyMRGCmosaic, mRGCMosaicSubspaceResponsesFileName);

        % Method to generate and save input cone mosaic subspace responses
        generateInputConeMosaicSubspaceRFmappingLinearResponses(...
            theRGCMosaic, opticsToEmploy, ...
            stimSizeDegs, stimXYpositionDegs, ...
            stimulusChromaticity, coneFundamentalsOptimizedForStimPosition, ...
            responsesFileName, varargin);

        % Method to compute subspace modulation responses of a cone mosaic under some optics
        [theConeMosaicSubspaceLinearModulationResponses, theConeMosaicNullResponses, ...
         HartleySpatialModulationPatterns, spatialSupportDegs, lIndices, mIndices] = ...
                   computeConeMosaicSubspaceRFmappingLinearResponses(theConeMosaic, theOptics,  ...
                                           thePresentationDisplay, ...
                                           stimParams, ...
                                           stimPositionDegs, ...
                                           varargin);

        % Method to generate and save input cone mosaic m-sequence responses
        generateInputConeMosaicMSequenceRFmappingLinearResponses(theRGCMosaic, opticsToEmploy, ...
            stimSizeDegs, stimXYpositionDegs, rfPixelsAcross, ...
            ternaryInsteadOfBinaryMsequence, mSequenceBitLength, ...
            stimulusChromaticity, coneFundamentalsOptimizedForStimPosition, ...
            responsesFileName, varargin);

        % Method to compute m-sequence modulation responses of a cone mosaic under some optics
        [theConeMosaicMSequenceLinearModulationResponses, theConeMosaicNullResponses, ...
         theConeMosaicMSequenceForwardModulationResponses, theConeMosaicMSequenceInverseModulationResponses, ...
         mSequenceIndicatorFunctions, spatialSupportDegs] = ...
                    computeConeMosaicMSequenceRFmappingLinearResponses(theConeMosaic, theOptics,  ...
                            thePresentationDisplay, ...
                            stimParams, ...
                            stimXYpositionDegs, ...
                            varargin);

        % Method to reset the parpool
        [shutdownParPoolOnceCompleted, numWorkers] = resetParPool(parPoolSize);

    end % Static methods

    % Main methods
    methods (Static)

        % Main method for doing mRGCmosaic analyses/computations
        run(restartPool);

        % Method to get the user's selection on which operation to perform
        operationSetToPerformContains = operationsMenu(mosaicParams);

        % Method to compute the compute the visually projected anatomical
        % cone Rcs
        performComputeVisuallyProjectedAnatomicalConeRcsByFittingRFmap(mosaicParams, varargin);

        % Method to perform the generateRGCMosaicOperation()
        performGenerateCenterConnectedRGCMosaicOp(mosaicParams, varargin);

        % Method to perform the
        % visualizeCenterConnectedRGCMosaicAndRemoveUnwantedRGCs operation
        performVisualizeCenterConnectedRGCMosaicAndRemoveUnwantedRGCsOp(mosaicParams, varargin);

        % Method to perform the visualizePSFsWithinRGCMosaic operation
        performVisualizePSFsWithinRGCMosaicOp(mosaicParams, tickSeparationArcMin);

        % Method to perform the visualizeConePSFsWithinRGCMosaic operation
        performVisualizeConePSFsWithinRGCMosaicOp(mosaicParams, tickSeparationArcMin);

        % Method to perform the computeInputConeMosaicSTFresponses operation
        performComputeInputConeMosaicSTFresponsesOp(mosaicParams);

        % Method to perform the computeInputConeMosaicVisuallyProjectedRcDegs operation
        performComputeInputConeMosaicVisuallyProjectedRcDegsOp(mosaicParams);

        % Method to perform the optimizeSurroundConePoolingModels operation
        performOptimizeSurroundConePoolingModelsOp(mosaicParams, varargin);

        % Method to perform the inspectOptimizedSurroundConePoolingModels operation
        performInspectOptimizedSurroundConePoolingModelsOp(mosaicParams, varargin);

        % Method to perform the summarizeOptimizedSurroundConePoolingModels operation
        performSummarizeOptimizedSurroundConePoolingModelsOp(mosaicEccsForSummaryStatistics, varargin);

        % Method to perform the generateComputeReadyMidgetRGCMosaic operation
        performGenerateComputeReadyMidgetRGCMosaicOp(mosaicParams);

        % Method to perform the computeVisualSTFsOfTheComputeReadyMidgetRGCmosaic operation
        performComputeVisualSTFsOfTheComputeReadyMidgetRGCMosaicOp(mosaicParams);

        % Method to perform the fitVisualSTFsOfTheComputeReadyMidgetRGCmosaic operation
        performFitVisualSTFsOfTheComputeReadyMidgetRGCMosaicOp(mosaicParams);

        % Method to perform the adjustGainOfComputeReadyMRGCMosaicBasedOnVisualSTFdata
        performAdjustGainOfComputeReadyMRGCMosaicBasedOnVisualSTFdata(mosaicParams);

        % Method to perform the restGainOfComputeReadyMRGCMosaic
        performResetGainOfComputeReadyMRGCMosaic(mosaicParams);

        % Method to perform the VisualizeConePoolingRFmapAndVisualSTFforTargetRGC operation
        performVisualizeConePoolingRFmapAndVisualSTFforTargetRGC(mosaicParams, varargin);

        % Method to perform the VisualizeDoGparamsOfVisualSTFsOfSingleMidgetRGCMosaic operation
        performVisualizeDoGparamsOfVisualSTFsOfSingleMidgetRGCMosaic(mosaicParams);

        % Method to perform the VisualizeDoGparamsOfVisualSTFsOfMultipleMidgetRGCMosaic operation
        performVisualizeDoGparamsOfVisualSTFsOfMultipleMidgetRGCMosaic(mosaicEccsToInclude);

        % Method to perform the ComputeVisualRFsUsingSubSpaceMapping operation
        performComputeVisualRFsUsingSubSpaceMapping(mosaicParams, varargin);

        % Method to perform the ComputeVisualRFsUsingMSequenceMapping operation
        performComputeVisualRFsUsingMSequenceMapping(mosaicParams, varargin);

        % Method to perform the VisualizeVisualRFmapForTargetRGC
        performVisualizeVisualRFmapForTargetRGC(mosaicParams, stimPositionDegs, rfMappingPixelMagnificationFactor, varargin);

        % Method to perform the ComputeVisualRFmapsOfRFcenterSubregions
        performComputeVisualRFcenterMapsViaDirectConvolutionWithPSF(varargin);

        % Method to perform the ContrastSTFsAcrossDifferentOpticsOrChromaticities operation
        performContrastSTFsAcrossDifferentChromaticities(...
            mosaicParams, rawFiguresRoot, scaledFiguresRoot, varargin);

        % Method to perform the ContrastMSequenceRFsAcrossDifferentOpticsOrChromaticities operation
        performContrastMSequenceRFsAcrossDifferentChromaticities(...
            mosaicParams, ...
            stimPositionDegs, ternaryInsteadOfBinaryMsequence, ...
            mSequenceBitLength, rfPixelsAcross,...
            rawFiguresRoot, scaledFiguresRoot, varargin);

        % Method to ask the user which mRGC mosaic to use for computing
        [mosaicHorizontalEccentricityDegs, mosaicEccsForSummaryStatistics] =  chooseMosaicToUse();
        
        % Method to query user and load desired computeReadyMRGCMosaic
        [theComputeReadyMRGCmosaic, fullPathToComputeReadyMRGCMosaic] = ...
            queryUserAndLoadComputeReadyMRGCMosaic(mosaicParams);

        % Method to ask the user which opticsParams to use for computing
        % the input cone mosaic STF responses
        [opticsParams, opticsToEmploy, inputConeMosaicSTFresponsesFileName] = ...
            chooseOpticsForInputConeMosaicSTFresponses(mosaicParams, varargin);

        % Method to ask the user what stimulus chromaticity to use for
        % computing the mosaicResponses and update the
        % mosaicResponsesFilename accordingly
        [stimulusChromaticity, coneFundamentalsOptimizedForStimPosition, mosaicResponsesFileName] = ...
            chooseStimulusChromaticityForMosaicResponsesAndUpdateFileName(...
                 mosaicResponsesFileName, identifierString, varargin);

        % Method to generate filenames for contrasted achromatic,
        % L- and M-cone isolating responses
        [mRGCMosaicAchromaticResponsesFileName, ...
          mRGCMosaicLconeIsolatingResponsesFileName, ...
          mRGCMosaicMconeIsolatingResponsesFileName, ...
          opticsParams, opticsToEmploy] = generateFileNamesForContrastedResponses(mosaicParams, measuredResponseType, userPrompt);

        % Method to encode m-sequence specific into to responses filename
        theFileName = encodeMSequenceSpecificInfoToFilename(theFileName, ...
            ternaryInsteadOfBinaryMsequence, mSequenceBitLength, rfPixelsAcross);

        % Method to ask the user which H1 cell index to use for optimizing
        % the RF surround cone pooling model
        [retinalRFmodelParams, gridSamplingScheme, optimizedRGCpoolingObjectsFileName] = ...
            chooseRFmodelForSurroundConePoolingOptimization(mosaicParams, opticsParams, varargin);

        % Method to specify the grid nodes for which we will optimize 
        % surround cone pooling 
        gridNodes = gridNodesToOptimize();

        % Method to generate filenames for various compute components
        [mosaicFileName, resourcesDirectory, pdfsDirectory] = ...
            resourceFileNameAndPath(component, varargin);

        % Method to load a pre-computed mRGC mosaic, solely based on its
        % hotizontal eccentricity
        theMRGCMosaic = loadPreComputedMRGCMosaic(horizontalEccDegs);

        % Method to obtain mosaicParams struct based on the mosaic's horizontal eccentricity
        mosaicParams = getMosaicParams(mosaicHorizontalEccentricityDegs);

        % Method to obtain the opticsParams based on the mosaicParams
        % Basically a wrapper for chooseOpticsForInputConeMosaicSTFresponses 
        % using the native optics and 0.0 refractive error
        opticsParams = getOpticsParams(mosaicParams);

        % Method to obtain the retinalRFmodelParams based on the
        % mosaicParams and opticsParams
        % Basically a wrapper for chooseRFmodelForSurroundConePoolingOptimization 
        % using the 4th H1 horizontal cell index
        retinalRFmodelParams = getSurroundParams(mosaicParams, opticsParams);

        % Various plotting methods
        plotRawCronerKaplanData();

        % M-sequence RF map smoothing ala ReidShapley (2002)
        [smoothedRF, smoothingKernel, rfPixelSizeSamples] = ...
            applyReidShapleySmoothingToRFmap(spatialSupportDegsX, theRFmap, rfPixelsAcross);

        % M-sequence RF map LUT ala ReidShapley (2002)
        [cLUToriginal, cLUTreversed] = generateReidShapleyRFmapLUT();

        % Scale figure and export it to a different file
        exportScaledFigure(sourcePDF, destinationPDF, scaleFactor);

    end

end
