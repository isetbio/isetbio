classdef MosaicPoolingOptimizer < handle

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

    end % Constants

    % Public methods (class interface)
    methods
        % Constructor
        function obj = MosaicPoolingOptimizer(theRGCmosaic, ...
            varargin)

            % Parse input
            p = inputParser;
            p.addRequired('theRGCmosaic', @(x)(isa(x, 'mRGCMosaic')));
            p.addParameter('minSpatialSamplingDegs', 0.25, @isnumeric);
            p.addParameter('samplingScheme', 'hexagonal', @(x)(ismember(x, {'hexagonal', 'rectangular'})));
            p.addParameter('visualizeSamplingGrids', false, @islogical);
            p.addParameter('generateSamplingGrids', false, @islogical);
            p.parse(theRGCmosaic, varargin{:});

            obj.theRGCMosaic = p.Results.theRGCmosaic;
            obj.retinalRFmodelParams = obj.defaultRetinalRFmodelParams;

            if (p.Results.generateSamplingGrids)
                fprintf('\nGenerating sampling grids. Please wait ...\n');
                % Generate the nominal multifocal spatial sampling grid
                obj.generateNominalSpatialSamplingGrid(p.Results.samplingScheme);
                
                % Generate the full multifocal sampling grids
                obj.generateSamplingGrids(p.Results.minSpatialSamplingDegs);
                fprintf('Done \n')

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

       
        % Method to load the computed input cone mosaic STF responses
        loadConeMosaicVisualSTFresponses(obj, responsesFileName);

        % Method to compute optimized RGC models (surround cone pooling weights)
        % for a grid node in the RGC mosaic
        compute(obj, gridNodeIndex, whichConeType, opticsParams, ...
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
            opticsParams, initialRetinalConePoolingParams, displayFittingProgress, figNo, figTitle);

        % Optimization components
        [modelConstants, retinalConePoolingParams, visualRcDegs] = computeOptimizationComponents(...
            obj, theRGCindex);
    end

    % Static methods
    methods (Static)

        % Method to generate the local dropbox path
        dropboxDir = localDropboxPath();

        % Fit a Gaussian model to a subregion STF
        [Gparams, theFittedSTF] = fitGaussianToSubregionSTF(...
            spatialFrequencySupportCPD, theSubregionSTF, ...
            RcDegsInitialEstimate, rangeForRc, multiStartsNum)

        % Fit a DoG model to an STF
        [DoGparams, theFittedSTF] = fitDifferenceOfGaussiansToSTF(...
             spatialFrequencySupportCPD, theSTF, ...
             RcDegsInitialEstimate, rangeForRc, multiStartsNum);

        % Method to select the highest-extending STF (across a set of STFs
        % measured at different orientations)
        theOptimalSTF = optimalSTFfromResponsesToAllOrientationsAndSpatialFrequencies(...
            orientationsTested, spatialFrequenciesTested, ...
            theResponsesAcrossAllOrientationsAndSpatialFrequencies)


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
            targetRGCposition, targetCenterConesNum, targetCenterConeMajorityType);

        % Method to visualize the visual RF map for a target RGC
        visualizeVisualRFmapForTargetRGC(...
            theComputeReadyMRGCmosaic, ...
            optimallyMappedSubspaceRFmapsFileName, ...
            targetRGCposition, targetCenterConesNum, targetCenterConeMajorityType, ...
            pdfFileName);

        % Method to visualize the visual RF maps for multiple target RGCs
        visualizeVisualRFmapsForMultipleTargetRGCs(...
            theComputeReadyMRGCmosaic, ...
            optimallyMappedSubspaceRFmapsFileName, ...
            mRGCMosaicSubspaceResponsesFileName, ...
            pdfFileName);

        % Method to setup the parameters and the display for conducting an
        % STF mapping experiment.
        [stimParams, thePresentationDisplay] = setupSTFmappingExperiment(inputConeMosaic, ...
            sceneFOVdegs, retinalImageResolutionDegs);


        % Method to compute visualSTF responses of a cone mosaic under
        % some optics
        [theConeMosaicSTFresponses, theConeMosaicNullResponses] = ...
            computeConeMosaicSTFresponses(theConeMosaic, theOptics, ...
                                       thePresentationDisplay, ...
                                       stimParams, stimPositionDegs, ...
                                       useParfor, visualizeResponses)

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

        % Method to visualize the fitted the visualSTF of the compute-ready MRGCmosaic
        visualizeFittedSTFsOfComputeReadyMidgetRGCMosaic(...
                computeReadyMosaicFilename, ...
                mRGCMosaicSTFresponsesFilename, ...
                employTemporalEquivalentEccentricity);

        
        visualizeConeMosaicSTFresponses(mRGCMosaicFileName, coneMosaicSTFresponsesFileName, varargin);


        % Method to setup the parameters and the display for conducting an
        % Subspace RF mapping experiment.
        [stimParams, thePresentationDisplay] = setupSubspaceRFmappingExperiment(wavelengthSupport, ...
             stimSizeDegs, maxSFcyclesPerDegree);

        % Method to compute visual RFs of the computeReadyMosaic (using the
        % subspace RF mapping method)
        computeVisualRFsOfComputeReadyMidgetRGCMosaic(...
            theComputeReadyMRGCmosaic, opticsParams, ...
            maxSFcyclesPerDegree, stimSizeDegs, posIncrementDegs, ...
            coneMosaicSubspaceResponsesFileName, ...
            mRGCMosaicSubspaceResponsesFileName, ...
            reComputeInputConeMosaicSubspaceRFmappingResponses, ...
            reComputeMRGCMosaicSubspaceRFmappingResponses, ...
            reComputeRFs, ...
            onlyVisualizeOptimallyMappedRFmaps, ...
            varargin);

        % Method to visualize the visual RF centers of all RGCs in a
        % compute ready MRGC mosaic
        visualizeVisualRFcentersOfComputeReadyMidgetRGCMosaic(...
            theComputeReadyMRGCmosaic, mRGCMosaicSubspaceResponsesFileName);

        % Method to generate and save input cone mosaic subspace responses
        generateInputConeMosaicSubspaceRFmappingLinearResponses(theRGCMosaic, theOptics, ...
            maxSFcyclesPerDegree, stimSizeDegs, stimXYpositionDegs, responsesFileName, varargin);

        % Method to compute subspace responses of a cone mosaic under some optics
        [theConeMosaicSubspaceLinearResponses, theConeMosaicNullResponses, ...
         HartleySpatialModulationPatterns, spatialSupportDegs, lIndices, mIndices] = ...
                   computeConeMosaicSubspaceRFmappingLinearResponses(theConeMosaic, theOptics,  ...
                                           thePresentationDisplay, ...
                                           stimParams, ...
                                           stimPositionDegs, ...
                                           varargin);

        % Method to reset the parpool
        [shutdownParPoolOnceCompleted, numWorkers] = resetParPool(parPoolSize);

    end % Static methods

    % Main methods
    methods (Static)

        % Main method for doing mRGCmosaic analyses/computations
        run();

        % Method to get the user's selection on which operation to perform
        operationSetToPerformContains = operationsMenu(mosaicParams);

        % Method to perform the generateRGCMosaicOperation()
        performGenerateCenterConnectedRGCMosaicOp(mosaicParams);

        % Method to perform the
        % visualizeCenterConnectedRGCMosaicAndRemoveUnwantedRGCs operation
        performVisualizeCenterConnectedRGCMosaicAndRemoveUnwantedRGCsOp(mosaicParams);

        % Method to perform the visualizePSFsWithinRGCMosaic operation
        performVisualizePSFsWithinRGCMosaicOp(mosaicParams);

        % Method to perform the computeInputConeMosaicSTFresponses operation
        performComputeInputConeMosaicSTFresponsesOp(mosaicParams);

        % Method to perform the optimizeSurroundConePoolingModels operation
        performOptimizeSurroundConePoolingModelsOp(mosaicParams);

        % Method to perform the inspectOptimizedSurroundConePoolingModels operation
        performInspectOptimizedSurroundConePoolingModelsOp(mosaicParams);

        % Method to perform the generateComputeReadyMidgetRGCMosaic operation
        performGenerateComputeReadyMidgetRGCMosaicOp(mosaicParams);

        % Method to perform the computeVisualSTFsOfTheComputeReadyMidgetRGCmosaic operation
        performComputeVisualSTFsOfTheComputeReadyMidgetRGCMosaicOp(mosaicParams);

        % Method to ask the user which opticsParams to use for computing
        % the inputConeMosaic STF responses
        [opticsParams, opticsToEmploy, coneMosaicSTFresponsesFileName] = ...
            chooseOpticsForInputConeMosaicSTFresponses(mosaicParams);

        % Method to ask the user which H1 cell index to use for optimizing
        % the RF surround cone pooling model
        [retinalRFmodelParams, gridSamplingScheme, optimizedRGCpoolingObjectsFileName] = ...
            chooseRFmodelForSurroundConePoolingOptimization(mosaicParams, opticsParams);

        % Method to specify the grid nodes for which we will optimize 
        % surround cone pooling 
        gridNodes = gridNodesToOptimize();

        % Method to generate filenames for various compute components
        [mosaicFileName, resourcesDirectory, pdfsDirectory] = ...
            resourceFileNameAndPath(component, varargin);
    end

end
