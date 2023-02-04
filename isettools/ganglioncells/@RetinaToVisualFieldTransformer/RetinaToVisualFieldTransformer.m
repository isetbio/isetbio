classdef RetinaToVisualFieldTransformer < handle

    % Public properties (read-only)
    properties (SetAccess = private)

       % The wavelength support for the employed PSF
       psfWavelengthSupport;

       % The @cMosaic object
       theConeMosaic;

       % Conversion factor from aperture diameter to characteristic radius
       coneCharacteristicRadiusConversionFactor;

       % The input optics params
       opticsParams;

       % The target visual params
       targetVisualRFDoGparams;

       % Whether to simulation the Croner&Kaplan estimation (i.e. fit line
       % weighting Gaussina functions to 1D projections of the visual 2D RF
       % maps.
       simulateCronerKaplanEstimation;

       % STF match mode: either match the shape or the DoG params of the
       % fitted shapes
       targetSTFmatchMode;

       % Whether to use a flatop gaussian model for estimating the characteristic
       % radius of the visual RF center
       flatTopGaussianForVisualRFcenterCharacteristicRadiusEstimation;

       % Subject ID for the chosen database and ranking
       testSubjectID;

       % Whether the subject requires subtraction of central refraction
       subtractCentralRefraction;

       % The spectrally-weighted PSFs
       theSpectrallyWeightedPSFData;

       % # of multi-starts
       multiStartsNum;
       multiStartsNumDoGFit;

       % Whether to do a dry run first
       doDryRunFirst;

       % The rfComputeStruct for L-center RFs
       LconeRFcomputeStruct;

       % The rfComputeStruct for M-center RFs
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
    end

    % Public methods (class interface)
    methods
        % Constructor
        function obj = RetinaToVisualFieldTransformer(theConeMosaic, ...
                opticsParams, targetVisualRFDoGparams, varargin) 
            % Parse input
            p = inputParser;
            p.addParameter('psfWavelengthSupport', [], @(x)(isvector(x)));
            p.addParameter('flatTopGaussianForVisualRFcenterCharacteristicRadiusEstimation', false, @islogical);
            p.addParameter('simulateCronerKaplanEstimation', true, @islogical);
            p.addParameter('targetSTFmatchMode', 'STFDoGparams', @(x)(ismember(x,{'STFDoGparams', 'STFcurve'})));
            p.addParameter('initialRetinalConePoolingParamsStruct', [], @(x)(isempty(x)||isstruct(x)));
            p.addParameter('multiStartsNumDoGFit', 64, @isscalar);
            p.addParameter('multiStartsNum', 10, @isscalar);
            p.addParameter('doDryRunFirst', false, @islogical);
            p.addParameter('computedRTVObjectExportDirectory', '', @(x)(isempty(x)||ischar(x)));
            p.parse(varargin{:});

            computedRTVObjectExportDirectory = p.Results.computedRTVObjectExportDirectory;
            initialRetinalConePoolingParamsStruct = p.Results.initialRetinalConePoolingParamsStruct;

            % Generate filename for saved object
            if (isempty(computedRTVObjectExportDirectory))
                fprintf('Computed object will NOT be saved to disk\n');
                obj.computedObjDataFileName = '';
            else
                obj.computedObjDataFileName = fullfile(computedRTVObjectExportDirectory, RetinaToVisualFieldTransformer.computedObjectDataFileName(opticsParams, targetVisualRFDoGparams));
                fprintf('Computed object will be saved to %s\n', obj.computedObjDataFileName);
            end

            obj.theConeMosaic = theConeMosaic;
            obj.coneCharacteristicRadiusConversionFactor = theConeMosaic.coneApertureToConeCharacteristicRadiusConversionFactor;

            obj.opticsParams = opticsParams;
            obj.targetVisualRFDoGparams = targetVisualRFDoGparams;

            obj.simulateCronerKaplanEstimation = p.Results.simulateCronerKaplanEstimation;
            obj.targetSTFmatchMode = p.Results.targetSTFmatchMode;
       
            obj.psfWavelengthSupport = p.Results.psfWavelengthSupport;
            obj.flatTopGaussianForVisualRFcenterCharacteristicRadiusEstimation = p.Results.flatTopGaussianForVisualRFcenterCharacteristicRadiusEstimation;
            obj.multiStartsNum = p.Results.multiStartsNum;
            obj.multiStartsNumDoGFit = p.Results.multiStartsNumDoGFit;
            obj.doDryRunFirst = p.Results.doDryRunFirst;

            % Assert that the cone mosaic contains the position
            % specified in the optical params struct
            coneMosaicOutline.x = theConeMosaic.eccentricityDegs(1) + theConeMosaic.sizeDegs(1)*0.5*[-1  1 1 -1 -1];
            coneMosaicOutline.y = theConeMosaic.eccentricityDegs(2) + theConeMosaic.sizeDegs(2)*0.5*[-1 -1 1  1 -1];
            in = inpolygon(obj.opticsParams.positionDegs(1), obj.opticsParams.positionDegs(2), ...
                coneMosaicOutline.x, coneMosaicOutline.y);
            assert(in, 'optics position (%2.3f,%2.3f) is outside the passed cone mosaic', ...
                obj.opticsParams.positionDegs(1), obj.opticsParams.positionDegs(2));

            % Assert that theConeMosaic and the optics refer to the same eye
            assert(strcmp(theConeMosaic.whichEye, obj.opticsParams.analyzedEye), ...
                'theConeMosaic and the optics params must having matching eye');

            % Assert that the centerConnectableConeTypes are valid cone types
            validConeTypes = [cMosaic.LCONE_ID cMosaic.MCONE_ID cMosaic.SCONE_ID];
            for iConeType = 1:numel(obj.targetVisualRFDoGparams.centerConnectableConeTypes)
                assert(ismember(obj.targetVisualRFDoGparams.centerConnectableConeTypes(iConeType), validConeTypes), ...
                       sprintf('Invalid centerConnectableConeType: %d.', obj.targetVisualRFDoGparams.centerConnectableConeTypes(iConeType)));
            end
            obj.targetVisualRFDoGparams.centerConnectableConeTypes = unique(obj.targetVisualRFDoGparams.centerConnectableConeTypes);

            % Assert that the surroundConnectableConeTypes are valid cone types
            validConeTypes = [cMosaic.LCONE_ID cMosaic.MCONE_ID cMosaic.SCONE_ID];
            for iConeType = 1:numel(obj.targetVisualRFDoGparams.surroundConnectableConeTypes)
                assert(ismember(obj.targetVisualRFDoGparams.surroundConnectableConeTypes(iConeType), validConeTypes), ...
                       sprintf('Invalid surroundConnectableConeType: %d.', obj.targetVisualRFDoGparams.surroundConnectableConeTypes(iConeType)));
            end
            obj.targetVisualRFDoGparams.surroundConnectableConeTypes = unique(obj.targetVisualRFDoGparams.surroundConnectableConeTypes);


            % Compute the spectrally-weighted PSFdata. This generates an
            % L-cone fundamental spectrally-weighted PSF, and an
            % M-cone fundamental spectrally-weighted PSF,
            % obj.theSpectrallyWeightedPSFData.LconeWeighted
            % obj.theSpectrallyWeightedPSFData.MconeWeighted
            % and
            % obj.theSpectrallyWeightedPSFData.LMconeWeighted
            obj.spectrallyWeightedPSFs();

            % Estimate the mean characteristic radius of cones at the examined position
            % as projected on to visual space using the computed PSF
            dStruct = obj.estimateConeCharacteristicRadiusInVisualSpace(...
                obj.opticsParams.positionDegs, ...
                obj.simulateCronerKaplanEstimation, ...
                targetVisualRFDoGparams.conesNumPooledByTheRFcenter);
    
            if (dStruct.conesNumInRetinalPatch==0)
                error('No cones in cone mosaic')
            end

            % Compute maxSpatialSupportDegs based on desired visual RF properties
            maxSpatialSupportDegs = dStruct.visualConeCharacteristicRadiusDegs * 6.0 * ...
                   sqrt(targetVisualRFDoGparams.conesNumPooledByTheRFcenter) * ...
                   targetVisualRFDoGparams.surroundToCenterRcRatio;
            maxSpatialSupportDegs = round(100*maxSpatialSupportDegs)/100;
            
            
            % Crop the PSFs to maxSpatialSupportDegs to speed up computations
            obj.cropPSF(maxSpatialSupportDegs);

            % Unit weights for center cones            
            indicesOfConesPooledByTheRFcenter = obj.targetVisualRFDoGparams.indicesOfConesPooledByTheRFcenter;  % dStruct.indicesOfConesSortedWithDistanceToTargetRFposition(1:targetVisualRFDoGparams.conesNumPooledByTheRFcenter);
            weightsOfConesPooledByTheRFcenter = obj.targetVisualRFDoGparams.weightsOfConesPooledByTheRFcenter;  

            % Compute retinal cone pooling params to generate the target
            % visual RF for L-cone center RGCs
            obj.LconeRFcomputeStruct = obj.retinalRFparamsForTargetVisualRF(indicesOfConesPooledByTheRFcenter, ...
                weightsOfConesPooledByTheRFcenter, targetVisualRFDoGparams, ...
                cMosaic.LCONE_ID, ...
                initialRetinalConePoolingParamsStruct.LconeRetinalConePoolingParams);

            % Compute retinal cone pooling params to generate the target
            % visual RF for M-cone center RGCs
            obj.MconeRFcomputeStruct = obj.retinalRFparamsForTargetVisualRF(indicesOfConesPooledByTheRFcenter, ...
                weightsOfConesPooledByTheRFcenter, targetVisualRFDoGparams, ...
                cMosaic.MCONE_ID, ...
                initialRetinalConePoolingParamsStruct.MconeRetinalConePoolingParams);

            % Save the computed object
            if (~isempty(computedRTVObjectExportDirectory))
                fprintf('Saving computed object to %s\n', obj.computedObjDataFileName);
                save(obj.computedObjDataFileName, 'obj');
            end

        end

        % Method to estimate the visually-projectected cone Rc given a target
        % position in the mosaic and corresponding PSF data
        dStruct = estimateConeCharacteristicRadiusInVisualSpace(obj, theTargetPositionDegs, ...
            simulateCronerKaplanEstimation, neighboringConesNum);

        % Method to remove large chunks of data
        % Called from midgetRGCMosaic.freeze()
        freeze(obj);

        % Method to overwrite the LconeRFcomputeStruct
        function overwriteLconeRFcomputeStruct(obj, newComputeStruct)
            obj.LconeRFcomputeStruct = newComputeStruct;
        end

        % Method to overwrite the MconeRFcomputeStruct
        function overwriteMconeRFcomputeStruct(obj, newComputeStruct)
            obj.MconeRFcomputeStruct = newComputeStruct;
        end

    end % Public methods

    % Private methods
    methods (Access=private)
        % Generate spectrally weighted PSFs
        spectrallyWeightedPSFs(obj);

        % Crop the PSFs
        cropPSF(obj,maxSpatialSupportDegs);

        % Method to compute the cone map for the RF center and its
        % corresponding Gaussian characteristic radius
        [visualRFcenterCharacteristicRadiusDegs, visualRFcenterConeMap, ...
         visualRFcenterCharacteristicRadiiDegs, visualRFcenterFlatTopExponents, ...
         visualRFcenterXYpos, visualRFcenterOrientationDegs, ...
         anatomicalRFcenterCharacteristicRadiusDegs] = analyzeRFcenter(obj, ...
            indicesOfConesPooledByTheRFcenter, weightsOfConesPooledByTheRFcenter, ...
            spatialSupportDegs, ...
            theEmployedPSFData);


        % Obtain theRFcomputeStruct (weights and indices of surround cones) by fitting the target visualRF
        theRFcomputeStruct = retinalRFparamsForTargetVisualRF(obj,...
            indicesOfConesPooledByTheRFcenter, ...
            weightsOfConesPooledByTheRFcenter, ...
            targetVisualRFDoGparams, ...
            centerConeType, ...
            initialRetinalConePoolingParamsStruct);
    
    end

    % Class methods
    methods (Static)
        % Method to generate a circularly symmetric  PSF from a given PSF
        theCircularPSF = circularlySymmetricPSF(thePSF, mode);

        % Center and rotate PSF
        data = centerAndRotatePSF(data);

       
        % Method to estimate various aspects of the geometry of a 2D shape
        [theCentroid, theAxesLengths, theRotationAngle] = estimateGeometry(supportX, supportY, zData);


        % Gaussian line weighting profile
        theLineWeightingProfile = gaussianLineWeightingProfile(params, spatialSupport);

        [visualConeCharacteristicRadiusDegs, bestHorizontalResolutionRotationDegs] = analyzeVisuallyProjectedConeAperture(theConeMosaic, ...
                 anatomicalConeApertureDiameterDegs, thePSFData, simulateCronerKaplanEstimation, hFig);


        % visual RF model: arbitrary shape (fixed) RF center and Gaussian surround
        [theRF, centerRF, surroundRF] = differenceOfArbitraryCenterAndGaussianSurroundRF(...
           modelConstants, paramsVector);

        % visual RF model: Gaussian center and Gaussian surround
        [theRF, centerRF, surroundRF] = differenceOfGaussianCenterAndGaussianSurroundRF(...
           modelConstants, paramsVector);

        % visual RF model: ellipsoidal Gaussian center and Gaussian surround
        [theRF, centerRF, surroundRF] = differenceOfEllipsoidalGaussianCenterAndGaussianSurroundRF(...
           modelConstants, paramsVector);

        % retinal cone pooling model: arbitrary center/double exponential surround
        pooledConeIndicesAndWeights = conePoolingCoefficientsForArbitraryCenterDoubleExpSurround(...
            modelConstants, conePoolingParamsVector);

        % retinal cone pooling model: arbitrary center/variable exponential surround
        pooledConeIndicesAndWeights = conePoolingCoefficientsForArbitraryCenterVariableExpSurround(...
            modelConstants, conePoolingParamsVector);

        % retinal cone pooling model: arbitrary center/gaussian surround
        pooledConeIndicesAndWeights = conePoolingCoefficientsForArbitraryCenterGaussianSurround(...
            modelConstants, conePoolingParamsVector);

        % retinal cone pooling model: arbitrary center/double gaussian surround
        pooledConeIndicesAndWeights = conePoolingCoefficientsForArbitraryCenterDoubleGaussianSurround(...
            modelConstants, conePoolingParamsVector);

        % retinal cone pooling model: arbitrary center/gaussian surround
        % with arbitrary adjustments in the surround
        pooledConeIndicesAndWeights = conePoolingCoefficientsForArbitraryCenterGaussianAdjustSurround(...
            modelConstants, conePoolingParamsVector);

        % Method to adjust the cone weights so as to counteract the gain
        % changes in relative cone efficacy
        [efficacyAdjustedConeWeights, maxEfficacy] = coneEfficacyAdjustedGains(theConeMosaic, ...
            coneDiametersDegs, outerSegmentLengthEccVariationAttenuationFactors, coneWeights, maxEfficacy);

        % Method to remove surround cone indices (and weights) for
        % non-connectable surround cones
        [surroundConeIndices, surroundConeWeights, ...
         nonConnectableSurroundConeIndices, nonConnectableSurroundConeWeights] = connectableSurroundConeIndicesAndWeights(...
            surroundConeIndices, surroundConeWeights, modelConstants)

        % Compute the fitted visualRF from the current retinal pooling params
        [theFittedVisualRF, theRetinalRFcenterConeMap, theRetinalRFsurroundConeMap, pooledConeIndicesAndWeights] = ...
            visualRFfromRetinalConePooling(modelConstants, retinalPoolingParams);

        % Method to compute the retinal subregion cone map by summing the
        % corresponding cone aperture maps
        retinalSubregionConeMap = retinalSubregionConeMapFromPooledConeInputs(theConeMosaic, ...
            conePosDegs, coneWeights, spatialSupportDegs);

        % Method to perform a simulation of the Crone&Kaplan STF measurement
        [theRMSEvector, theRotatedRF, theRFprofile, ...
          theVisualSTF, theSpatialFrequencySupport, ...
          theFittedSTFsurroundToCenterRcRatio, ...
          theFittedSTFsurroundToCenterIntegratedSensitivityRatio, ratioWeights] = performCronerKaplanSimulation(...
            theVisualRF, theRetinalRFcenterConeMap, theRetinalRFsurroundConeMap, ...
            currentRetinalPoolingParams, retinalConePoolingParams, ... 
            bestHorizontalResolutionRotationDegs, visualRFcenterCharacteristicRadiusDegs, ...
            sfSupport, targetSTFmatchMode, theTargetSTF, ...
            targetRsRcRatio, targetIntSensSCRatio, ...
            multiStartsNumDoGFit, figNo, figureName);

        % Method to fit a 2D Gaussian ellipsoid to a continuous RF map
        theFittedGaussian = fitGaussianEllipsoid(supportX, supportY, theRF, varargin);

         % Method to fit a 2D Gaussian ellipsoid to a scattered RF map
        theFittedGaussian = fitScatterGaussianEllipsoid(supportX, supportY, theRF, inputWeights, inputPositions, varargin);

        % Method to fit a Gaussian line weight function to a 1D RFmap profile
        theFittedGaussianLineWeightingFunction = fitGaussianLineWeightingFunction(...
            supportXdegs, theRFprofile, varargin);

        % Method to fit a Difference of Gaussians model to a spatial
        % transfer function
        [DoGparams, theFittedSTF] = fitDoGmodelToMeasuredSTF(...
            spatialFrequencySupport, theSTF, RcDegsInitialEstimate, multiStartsNum, varargin);

        % Method to detemine best RF rotation to maximize resolution along
        % horizontal axis, and rotate the RF according to this angle
        [theRotatedRF, bestHorizontalResolutionRotationDegs] = ...
            bestHorizontalResolutionRFmap(theRF, bestHorizontalResolutionRotationDegs);

        % Method to compute the one-sided STF from the RF profile
        [oneSidedSpatialFrequencySupport, oneSidedSTF] = spatialTransferFunction(spatialSupportDegs, theRFprofile);

        %Given Kc (center gain), KsToKcRatio (surround/center peak sensitivity
        % ratio) and Rwide (wide RF surround radius), narrowToWideVolumeRatio and RnarrowToRwideRatio,
        % compute the peak sensitivity gains for the wide and the narrow surround components
        [Kwide, Knarrow, Rnarrow] = H1doubleExponentRFparams(Kc, Rwide, KsToKcPeakRatio, narrowToWideVolumeRatio, RnarrowToRwideRatio);

        % Method to test the doule exponent RF params
        testDoubleExponentRFparams();

        % Method to visualize the fitted param values
        visualizeRetinalSurroundModelParametersAndRanges(ax, modelParams);

        % Method to generate the datafilename where the computed object is saved
        dataFileName = computedObjectDataFileName(opticsParams, targetVisualDoGparams);
    end


end
