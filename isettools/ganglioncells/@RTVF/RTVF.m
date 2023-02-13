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

       % The computed L-cone RF compute struct
       LconeRFcomputeStruct;

       % The computed M-cone RF compute struct
       MconeRFcomputeStruct;

       % The data filename where the computed object is saved
       computedObjDataFileName;
    end

    % Constant properties
    properties (Constant, Hidden)
        Artal = 'Artal2012';
        Polans = 'Polans2015';
    end

    
    properties (Constant)
        validZernikeDataBases = {...
            RetinaToVisualFieldTransformer.Artal ...
            RetinaToVisualFieldTransformer.Polans ...
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
                varargin) 

            p = inputParser;
            p.addParameter('computedRTVFobjectExportDirectory', '', @(x)(isempty(x)||ischar(x)));
            p.addParameter('visualizeSpectrallyWeightedPSFs', false, @islogical);
            p.addParameter('computeLconeCenterComputeStruct', true, @islogical);
            p.addParameter('computeMconeCenterComputeStruct', true, @islogical);
            p.addParameter('initialRetinalConePoolingParamsStruct', [], @(x)(isempty(x)||isstruct(x)));
            p.addParameter('multiStartsNumRetinalPooling', 1, @isscalar);
            p.addParameter('multiStartsNumDoGFit', 64, @isscalar);
            p.parse(varargin{:});

            % Handle optional inputs
            obj.multiStartsNumRetinalPooling = p.Results.multiStartsNumRetinalPooling;
            obj.multiStartsNumDoGFit = p.Results.multiStartsNumDoGFit;

            visualizeSpectrallyWeightedPSFs = p.Results.visualizeSpectrallyWeightedPSFs;
            initialRetinalConePoolingParamsStruct = p.Results.initialRetinalConePoolingParamsStruct;
            
            computedRTVFobjectExportDirectory = p.Results.computedRTVFobjectExportDirectory;
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
                fprintf('Computed object will NOT be saved to disk.\n');
                obj.computedObjDataFileName = '';
            else
                obj.computedObjDataFileName = fullfile(computedRTVFobjectExportDirectory, obj.computeObjectDataFileName());
                fprintf('Computed object will be saved to %s\n', obj.computedObjDataFileName);
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
            maxRsRcRatio =  max([ ...
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
%             if (computeLconeCenterComputeStruct && computeMconeCenterComputeStruct)
%                 % Do them in parallel
%                 computeStructs = cell(1,2)
%                 parfor i = 1:2
%                 end
% 
%             else
                if (computeLconeCenterComputeStruct)
                    obj.LconeRFcomputeStruct = obj.retinalConePoolingParamsForTargetVisualRF(...
                        cMosaic.LCONE_ID, ...
                        initialLconeRetinalConePoolingParams);
                end
    
                if (computeMconeCenterComputeStruct)
                    obj.MconeRFcomputeStruct = obj.retinalConePoolingParamsForTargetVisualRF(...
                        cMosaic.MCONE_ID, ...
                        initialMconeRetinalConePoolingParams);
                end
            %end


            % Save the computed object
            obj.saveComputedObject(computeLconeCenterComputeStruct, computeMconeCenterComputeStruct);

        end % 
        % Constructor
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

        % Method to compute the visual RF ensuing from the current retinal
        % cone pooling params
        [theVisualRF, theRetinalRFcenterConeMap, theRetinalRFsurroundConeMap, pooledConeIndicesAndWeights] = ...
            visualRFfromRetinalConePooling(obj, modelConstants, retinalPoolingParams)

        % Method to compute the Croner&Kaplan RF analysis
        dataOut = visualRFmapPropertiesFromCronerKaplanAnalysis(obj, theVisualRF);

        % Method to compute the RFcomputeStruct
        RFcomputeStruct = retinalConePoolingParamsForTargetVisualRF(obj, ...
                centerConeType, initialRetinalConePoolingParamsStruct);

        % Method to compute the saved obj filename
        dataFileName = computeObjectDataFileName(obj);

        % Method to save the computed RTVF object (or parts of it)
        saveComputedObject(obj,computeLconeCenterComputeStruct, computeMconeCenterComputeStruct);
    end % private methods


    % Class methods
    methods (Static)
        % Method to determine the centroid,and minor/major axis lengehs of a map
        %[theCentroid, theAxesLengths] = estimateGeometry(supportX, supportY, rfMap)

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
        displayFittedModel(figNo, modelConstants, theVisualRFmap, theVisualSTFdata, ...
    targetParams);

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

        % Method to inspect a saved RTVF object
        generateSpatialRFs(theRTVFobj);
    end

end

