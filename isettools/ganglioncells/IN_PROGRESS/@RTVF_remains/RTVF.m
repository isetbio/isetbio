classdef RTVF < handle

    % Public properties (read-only)
    properties (SetAccess = private)
       % The @cMosaic object
       coneMosaic;

       % The input optics params
       opticsParams;

       % The target visual params
       targetVisualRFDoGparams;

       % L- and M-cone weighted PSFs 
       % appropriate for the optical
       % position 
       spectrallyWeightedPSFData;

       % Cone aperture blur kernel
       coneApertureBlurKernel;

       % The visual RF center of the retinal cone map
       visualRFcenterConeMap;

       % The rotation that results in the narrowest x-profile for the
       % visual RF center
       bestHorizontalResolutionRotationDegs;

       % The characteristic radius of the visual RF center as estimated by
       % fitting a 1D Gaussian line weighting function to the 1D profile
       % of the visual RF center
       visualRFcenterRcDegs;

       % The characteristic radius of the anatomical RF center
       anatomicalRFcenterCharacteristicRadiusDegs;

       % How many multistarts to use for fitting a DoG model to the
       % visualRF obtained with the current retinal cone pooling params
       multiStartsNumDoGFit;

       % How many multistarts to use for the retinal cone pooling params
       multiStartsNumRetinalPooling;
       useParallelMultiStart;

       % Compute method
       stfComputeMethod;

       % Compute method resources
       stfComputeMethodResources;

       % The computed L-cone RF compute struct
       LconeRFcomputeStruct;

       % The computed M-cone RF compute struct
       MconeRFcomputeStruct;

       % The data filename where the computed object is saved
       computedObjDataFileName;

       % Boolean indicating whether the RTVF file was updates
       RTVFfileWasUpdated;
    end

    % Constant properties
    properties (Constant, Hidden)
        Artal = 'Artal2012';
        Polans = 'Polans2015';

        directSTFcomputeMethod = 'cone mosaic STF response based';
        modeledSTFcomputeMethod = 'modeled STF';
    end

    
    properties (Constant)
        validSTFcomputeMethods = {
            RTVF.directSTFcomputeMethod ...
            RTVF.modeledSTFcomputeMethod ...
            };

        validZernikeDataBases = {...
            RTVF.Artal ...
            RTVF.Polans ...
            };

        ArtalSubjectsNum = 54;
        PolansSubjectsNum = 10;

        defaultOpticsParams = struct(...
            'ZernikeDataBase', 'Polans2015', ...
            'examinedSubjectRankOrder', 6, ...
            'pupilDiameterMM', 3.0, ...
            'positionDegs', [0 0], ...
            'refractiveErrorDiopters', 0.0, ...    % use -999 for optics that do not subtract the central refraction
            'analyzedEye', 'right eye', ...
            'subjectRankingEye', 'right eye', ...
            'psfUpsampleFactor', 1.0, ...
            'wavefrontSpatialSamples', 401);

        % Default defaultTargetVisualRFDoGParams: no surround
        defaultTargetVisualRFDoGParams = struct(...
            'retinalConePoolingModel', 'arbitraryCenterConeWeights_doubleExpH1cellIndex1SurroundWeights', ...
            'centerConnectableConeTypes', [cMosaic.LCONE_ID cMosaic.MCONE_ID], ...
            'surroundConnectableConeTypes', [cMosaic.LCONE_ID cMosaic.MCONE_ID], ...
            'indicesOfConesPooledByTheRFcenter', [], ...
            'weightsOfConesPooledByTheRFcenter', [], ...
            'surroundToCenterRcRatio', 0.0, ...
            'surroundToCenterIntegratedSensitivityRatio', 0);
    end

    % Public methods (class interface)
    methods
        % Constructor
        function obj = RTVF(...
                theConeMosaic, ...
                theOpticsParams, ...
                theTargetVisualRFDoGparams, ...
                stfComputeMethod, ...
                varargin) 

            p = inputParser;
            p.addParameter('computedRTVFobjectExportDirectory', '', @(x)(isempty(x)||ischar(x)));
            p.addParameter('progressFigureName', '', @ischar);
            p.addParameter('visualizeSpectrallyWeightedPSFs', false, @islogical);
            p.addParameter('computeMethodResources', @isstruct);
            p.addParameter('computeLconeCenterComputeStruct', true, @islogical);
            p.addParameter('computeMconeCenterComputeStruct', true, @islogical);
            p.addParameter('initialRetinalConePoolingParamsStruct', [], @(x)(isempty(x)||isstruct(x)));
            p.addParameter('multiStartsNumRetinalPooling', 1, @isscalar);
            p.addParameter('useParallelMultiStart', false, @islogical);
            p.addParameter('multiStartsNumDoGFit', 64, @isscalar);
            p.parse(varargin{:});

            % Handle optional inputs
            obj.multiStartsNumRetinalPooling = p.Results.multiStartsNumRetinalPooling;
            obj.multiStartsNumDoGFit = p.Results.multiStartsNumDoGFit;
            obj.useParallelMultiStart = p.Results.useParallelMultiStart;
            
            visualizeSpectrallyWeightedPSFs = p.Results.visualizeSpectrallyWeightedPSFs;
            initialRetinalConePoolingParamsStruct = p.Results.initialRetinalConePoolingParamsStruct;
            
            computedRTVFobjectExportDirectory = p.Results.computedRTVFobjectExportDirectory;
            progressFigureName = p.Results.progressFigureName;
            
            obj.stfComputeMethod = stfComputeMethod;
            obj.stfComputeMethodResources = p.Results.computeMethodResources;

            computeLconeCenterComputeStruct = p.Results.computeLconeCenterComputeStruct;
            computeMconeCenterComputeStruct = p.Results.computeMconeCenterComputeStruct;

            % Handle empty opticsParams
            if (isempty(theOpticsParams))
                obj.opticsParams = RTVF.defaultOpticsParams;
            else
                obj.opticsParams = theOpticsParams;
            end

            % Handle empty cone mosaic
            if (isempty(theConeMosaic))
                obj.generateDefaultConeMosaic(theOpticsParams.positionDegs);
            else
                obj.coneMosaic = theConeMosaic;
            end
            obj.coneMosaic.visualize();

            % Handle empty argetVisualRFDoGparams
            if (isempty(theTargetVisualRFDoGparams))
                obj.targetVisualRFDoGparams = RTVF.defaultTargetVisualRFDoGParams;

                % Median value of C&K Rs/Rc ratios
                [temporalEccDegs, RcRsRatios] = ...
                    RGCmodels.CronerKaplan.digitizedData.parvoCenterSurroundRadiusRatioAgainstEccentricity();
                obj.targetVisualRFDoGparams.surroundToCenterRcRatio = prctile(1./RcRsRatios, 50);
            

                % Median value of C&K S/C ratios
                [temporalEccDegs, SCratios] = ...
                    RGCmodels.CronerKaplan.digitizedData.parvoSurroundCenterIntSensisitivityRatioAgainstEccentricity();
                obj.targetVisualRFDoGparams.surroundToCenterIntegratedSensitivityRatio = prctile(SCratios,50);
            else
                obj.targetVisualRFDoGparams = theTargetVisualRFDoGparams;
            end
            

            % Generate filename for saved object
            if (isempty(computedRTVFobjectExportDirectory))
                fprintf('\nComputed object will NOT be saved to disk.\n');
                obj.computedObjDataFileName = '';
            else
                obj.computedObjDataFileName = fullfile(computedRTVFobjectExportDirectory, obj.computeObjectDataFileName());
                fprintf('\nComputed object will be saved to %s\n', obj.computedObjDataFileName);
            end

            % Compute spectrally weighted (L-cone, M-cone and L+M-cone weighted)
            %         PSFs, where the L- and M-w spectral eights are derived from the retinal
            %         spectral absorptance of actual L- and M-cones in the input cone
            %         mosaic at the desrired eccentricity (opticsParams.positionDegs), 
            %         taking into account macular pigment density at that location, as 
            %         specified by the input @cMosaic.
            %         This sets the obj.theSpectrallyWeightedPSFData parameter
            %         and updates the obj.opticsParams
            obj.computeSpectrallyWeightedPSFs(visualizeSpectrallyWeightedPSFs);

            % Generate the cone aperture blur kernel for the RF center pooled cones
            %          i.e., targetVisualRFDoGparams.indicesOfConesPooledByTheRFcenter
            %          This sets the obj.coneApertureKernel parameter and
            %          it is used to compute retinal cone maps for both the
            %          center and the surround subregions
            obj.computeConeApertureBlurKernel();

            % Compute the characteristic radius of the retinal RFcenter cone map
            %          as projected in visual space using the computed L+M-cone weighted PSF.
            [obj.anatomicalRFcenterCharacteristicRadiusDegs, ...
             obj.visualRFcenterRcDegs, ...
             obj.visualRFcenterConeMap, ...
             obj.bestHorizontalResolutionRotationDegs] = obj.rfCenterCharacteristicRadiusInVisualSpace();


            % Retrieve the C&K Rs/Rc ratios
            [CronerKaplanRcRsRatios.temporalEccDegs, CronerKaplanRcRsRatios.val] = ...
                RGCmodels.CronerKaplan.digitizedData.parvoCenterSurroundRadiusRatioAgainstEccentricity();

            % Max Rs/Rc ratio (80% prctile of the C&K data)
            maxRsRcRatio = max([ ...
                prctile(1./CronerKaplanRcRsRatios.val, 80), ...
                obj.targetVisualRFDoGparams.surroundToCenterRcRatio]);

            % Crop the PSF to this max spatial support
            obj.cropPSF(obj.visualRFcenterRcDegs * maxRsRcRatio);

            % Recompute the visual RFcenter cone map which may have
            % different spatial support based on the above cropping
            [~, ~, obj.visualRFcenterConeMap, obj.bestHorizontalResolutionRotationDegs] = ...
             obj.rfCenterCharacteristicRadiusInVisualSpace( ...
                    'computeAnatomicalRFcenterRc', false, ...
                    'computeVisualRFcenterRc', false);

            % Initial retinal cone pooling params
            if (~isempty(initialRetinalConePoolingParamsStruct))
                if (isfield(initialRetinalConePoolingParamsStruct, 'LconeRFcenter'))
                    initialLconeRetinalConePoolingParams = initialRetinalConePoolingParamsStruct.LconeRFcenter;
                end
                if (isfield(initialRetinalConePoolingParamsStruct, 'MconeRFcenter'))
                    initialMconeRetinalConePoolingParams = initialRetinalConePoolingParamsStruct.MconeRFcenter;  
                end
            else
                initialLconeRetinalConePoolingParams = [];
                initialMconeRetinalConePoolingParams = [];
            end

            % Action!
            if (strcmp(obj.stfComputeMethod, RTVF.modeledSTFcomputeMethod))
                if (computeLconeCenterComputeStruct) 
                    obj.LconeRFcomputeStruct = obj.retinalConePoolingParamsForTargetVisualRF(...
                        cMosaic.LCONE_ID, ...
                        initialLconeRetinalConePoolingParams, ...
                        progressFigureName);
                end
    
                if (computeMconeCenterComputeStruct)
                    obj.MconeRFcomputeStruct = obj.retinalConePoolingParamsForTargetVisualRF(...
                        cMosaic.MCONE_ID, ...
                        initialMconeRetinalConePoolingParams, ...
                        progressFigureName);
                end
            else
                % The first majority cone type
                switch (obj.targetVisualRFDoGparams.targetRGCmajorityConeType)
                    case cMosaic.LCONE_ID
                        initialRetinalConePoolingParams = initialLconeRetinalConePoolingParams;
                        fprintf(2,'Computing L-cone compute struct\n');
                    case cMosaic.MCONE_ID
                        initialRetinalConePoolingParams= initialMconeRetinalConePoolingParams;
                        fprintf(2,'Computing M-cone compute struct\n');
                    otherwise
                        initialRetinalConePoolingParams = initialLconeRetinalConePoolingParams;
                end

                theComputeStruct = obj.retinalConePoolingParamsForTargetVisualRF(...
                        [], ...
                        initialRetinalConePoolingParams, ...
                        progressFigureName);

                switch (obj.targetVisualRFDoGparams.targetRGCmajorityConeType)
                    case cMosaic.LCONE_ID
                        obj.LconeRFcomputeStruct = theComputeStruct;
                    case cMosaic.MCONE_ID
                        obj.MconeRFcomputeStruct = theComputeStruct;
                    otherwise
                        obj.LconeRFcomputeStruct = theComputeStruct;
                        obj.MconeRFcomputeStruct = theComputeStruct;
                end

                
                % Now the other majory cone type
                obj.targetVisualRFDoGparams.indicesOfConesPooledByTheRFcenter = ...
                    obj.targetVisualRFDoGparams.indicesOfConesPooledByTheRFcenterOfDifferentMajorityConeType;
                obj.targetVisualRFDoGparams.weightsOfConesPooledByTheRFcenter = ...
                    obj.targetVisualRFDoGparams.weightsOfConesPooledByTheRFcenterOfDifferentMajorityConeType;

                switch (obj.targetVisualRFDoGparams.targetRGCdifferentMajorityConeType)
                    case cMosaic.LCONE_ID
                        initialRetinalConePoolingParams = initialLconeRetinalConePoolingParams;
                        fprintf(2,'Computing L-cone compute struct\n');
                    case cMosaic.MCONE_ID
                        initialRetinalConePoolingParams = initialMconeRetinalConePoolingParams;
                        fprintf(2,'Computing M-cone compute struct\n');
                    otherwise
                        initialRetinalConePoolingParams = initialLconeRetinalConePoolingParams;
                end

                theComputeStruct = obj.retinalConePoolingParamsForTargetVisualRF(...
                        [], ...
                        initialRetinalConePoolingParams, ...
                        progressFigureName);

                switch (obj.targetVisualRFDoGparams.targetRGCdifferentMajorityConeType)
                    case cMosaic.LCONE_ID
                        obj.LconeRFcomputeStruct = theComputeStruct;
                    case cMosaic.MCONE_ID
                        obj.MconeRFcomputeStruct = theComputeStruct;
                    otherwise
                        obj.LconeRFcomputeStruct = theComputeStruct;
                        obj.MconeRFcomputeStruct = theComputeStruct;
                end
            end

            % Remove the stfComputeMethodResources. No longer needed after the
            % computation is done
            obj.removeResources();

            % Save the computed object
            if (isempty(computedRTVFobjectExportDirectory))
                obj.RTVFfileWasUpdated = false;
                fprintf('Computed RTVF object was not saved to the disk.');
            else
                obj.RTVFfileWasUpdated  = obj.saveComputedObject(computeLconeCenterComputeStruct, computeMconeCenterComputeStruct);
            end
        end % 
        % Constructor

        function removeResources(obj)
            if (strcmp(obj.stfComputeMethod, RTVF.directSTFcomputeMethod))
                obj.stfComputeMethodResources = [];
            end
        end

        % Method to compute the visual RF ensuing from the current retinal
        % cone pooling params
        [theVisualRF, theRetinalRFcenterConeMap, theRetinalRFsurroundConeMap, pooledConeIndicesAndWeights] = ...
            visualRFfromRetinalConePooling(obj, modelConstants, retinalPoolingParams)

        % Method to compute the visual STF from the visualRF simulating the Croner&Kaplan RF analysis
        dataOut = visualSTFfromCronerKaplanAnalysisOfVisualRF(obj, theVisualRF, recomputeBestHorizontalResolutionRFmap);

<<<<<<< HEAD:isettools/ganglioncells/IN_PROGRESS/@RTVF_remains/RTVF.m
        % Method to compute the visual STF from cone mosaic STF responses simulating the Croner&Kaplan RF analysis
        dataOut = visualSTFfromCronerKaplanAnalysisOfconeMosaicSTFresponses(obj, pooledConeIndicesAndWeights);

=======
>>>>>>> master:isettools/ganglioncells/@RTVF/RTVF.m
        % Method to freeze the obj (i.e., remove large chunks of data that
        % are of no use after the cone weights to the surround have been
        % computed)
        freeze(obj);
    end % public methods


    % Private methods
    methods (Access=private)

        % PSF computation
        computeSpectrallyWeightedPSFs(obj, visualize, varargin);
        cropPSF(obj, maxSpatialSupportDegs);

        % Subregion computation
        computeConeApertureBlurKernel(obj);
        retinalSubregionConeMap = retinalSubregionMapFromPooledConeInputs(obj, ...
            conePosDegs, coneApertureAreas, osLengthAttenuationFactors, coneWeights);

        % RF center params
        [anatomicalRFcenterCharacteristicRadiusDegs, ...
         visualRFcenterRcDegs, ...
         visualRFcenterConeMap, ...
         bestHorizontalResolutionRotationDegs] = rfCenterCharacteristicRadiusInVisualSpace(obj, varargin);

        % Method to compute the RFcomputeStruct
        RFcomputeStruct = retinalConePoolingParamsForTargetVisualRF(obj, ...
                centerConeType, initialRetinalConePoolingParamsStruct, ...
                progressFigureName);

        % Method to compute the saved obj filename
        dataFileName = computeObjectDataFileName(obj);

        % Method to save the computed RTVF object (or parts of it)
        RTVFfileWasUpdated  = saveComputedObject(obj,computeLconeCenterComputeStruct, computeMconeCenterComputeStruct);
    end % private methods


    % Class methods
    methods (Static)
        % Method to compute the retinal spectral quantal efficiencies 
        % at a target position within a @cMosaic
        % taking into account variations in MP density
        [L,M,wavelengthSupport] = LMconeSpectralWeightings(theConeMosaic, theTargetEccDegs)

        % Method to determine the centroid,and minor/major axis lengehs of a map
        [theCentroid, theAxesLengths, theRotationAngle] = estimateGeometry(supportX, supportY, rfMap)

        % Method to find the rotation of theRF which results in the best
        % resolution (narrowest profile) along the x-axis
        [theRotatedRF, rotationDegs] = bestHorizontalResolutionRFmap(theRF, rotationDegs, debugRadonTransformAnalysis);

        % Method to fit a 2D Gaussian ellipsoid to a RF cone map
        theFittedGaussian = fitGaussianEllipsoid(supportX, supportY, theRFconeMap, varargin);

        % Method to fit the RF line-weighting profile with a Gaussian
        theFittedGaussianLineWeightingFunction = ...
            fitGaussianLineWeightingFunction(spatialSupport, theRFprofile, varargin)

        % Method to compute the line weighting function for a Gaussian
        theLineWeightingProfile = gaussianLineWeightingProfile(params, spatialSupport);

        % Method to compute the STF  of a RF line weighting function
        [oneSidedSpatialFrequencySupport, oneSidedSTF] = spatialTransferFunction(...
            spatialSupportDegs, theLineWeightingFunction);

        % Method to fit a DoG STF to the measured STF
        [DoGparams, theFittedSTF] = fitDoGmodelToMeasuredSTF(...
            sfCPD, theMeasuredSTF, RcDegsInitialEstimate, rangeForRc, multiStartsNum);

        % Method to visualize the params and ranges of a fitted model
        visualizeFittedModelParametersAndRanges(ax, modelParams);

        % Method to visualize the fitted model
        displayFittedModel(figNo, modelConstants, theVisualRFmap, theVisualSTFdata, targetParams);

        % Method to visualize the fitting progress
        rmseSequence = displayFittingProgress(hFigProgress, videoOBJ, rmseSequence, ...
                RsRcRatioResidual, SCintSensRatioResidual, theCurrentRMSE, ...
                retinalConePoolingParams, currentRetinalPoolingParamValues, ...
                theCurrentSTFdata, ...
                spatialSupportDegsX, ...
                spatialSupportDegsY, ...
                theCurrentRetinalRFcenterConeMap, ...
                theCurrentRetinalRFsurroundConeMap);

        % Method to compute the pooling cone indices and weights for the
        % current conePoolingParamsVector for the
        % 'arbitraryCenterConeWeights_doubleExpH1cellIndex1SurroundWeights'
        % retinal cone pooling model
        pooledConeIndicesAndWeights = conePoolingCoefficientsForArbitraryCenterDoubleExpSurround( ...
            modelConstants, conePoolingParamsVector)

        % Method to compute H1 double exponent RF params that are dependent
        % on other H1 model params
        [Kwide, Knarrow, Rnarrow] = H1doubleExponentRFparams(...
            Kc, Rwide, KsToKcPeakRatio, narrowToWideVolumeRatio, RnarrowToRwideRatio);

        % Method to compute the surround cone indices and weights based on
        % the the obj.targetVisualRFDoGparams.surroundConnectableConeTypes
        % (which is encoded in the modelConstant struct)
        [surroundConeIndices, surroundConeWeights, ...
          nonConnectableSurroundConeIndices, ...
          nonConnectableSurroundConeWeights] = connectableSurroundConeIndicesAndWeights(...
                surroundConeIndices, surroundConeWeights, modelConstants);

        % Method to generate the center/surround spatial RF from an RTVFobj
        generateSpatialRFs(theRTVFobj);

        % Method to accept either a 'y' or a 'n' response
        val = queryUserForYesNoResponse(message);
    end

end

