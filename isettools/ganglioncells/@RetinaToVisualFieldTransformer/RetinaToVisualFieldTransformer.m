classdef RetinaToVisualFieldTransformer < handle

    % Public properties (read-only)
    properties (SetAccess = private)

       % The wavelength support for the employed PSF
       psfWavelengthSupport;

       % The PSF circular symmetry mode
       psfCircularSymmetryMode;

       % The @cMosaic object
       theConeMosaic;

       % Conversion factor from aperture diameter to characteristic radius
       coneCharacteristicRadiusConversionFactor;

       % The input optics params
       opticsParams;

       % The target visual params
       targetVisualRFDoGparams;

       % Whether to use a flatop gaussian model for estimating the characteristic
       % radius of the visual RF center
       flatTopGaussianForVisualRFcenterCharacteristicRadiusEstimation;

       % Subject ID for the chosen database and ranking
       testSubjectID;

       % Whether the subject requires subtraction of central refraction
       subtractCentralRefraction;

       % The computed (vLambda weighted, or not) PSFs
       thePSFData;
       theCircularPSFData;

       % # of multi-starts
       multiStartsNum;

       % Whether to do a dry run first
       doDryRunFirst;

       % The results of the computation
       rfComputeStruct;

       % The data filename where the computed object is saved
       computedObjDataFileName;
    end

    % Constant properties
    properties (Constant, Hidden)
        Artal = 'Artal2012';
        Polans = 'Polans2015';
        psfCircularSymmetryModeNone = 'none';
        psfCircularSymmetryModeAverage = 'average';
        psfCircularSymmetryModeBestResolution = 'bestResolution';
        psfCircularSymmetryModeWorstResolution = 'worstResolution';
    end

    properties (Constant)
        validZernikeDataBases = {...
            RetinaToVisualFieldTransformer.Artal ...
            RetinaToVisualFieldTransformer.Polans ...
            };
        validPSFcircularSymmetryModes = {...
            RetinaToVisualFieldTransformer.psfCircularSymmetryModeNone ...
            RetinaToVisualFieldTransformer.psfCircularSymmetryModeAverage ...
            RetinaToVisualFieldTransformer.psfCircularSymmetryModeBestResolution ...
            RetinaToVisualFieldTransformer.psfCircularSymmetryModeWorstResolution ...
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
            p.addParameter('psfCircularSymmetryMode', ...
                RetinaToVisualFieldTransformer.psfCircularSymmetryModeBestResolution, ...
                @(x)(ismember(x, RetinaToVisualFieldTransformer.validPSFcircularSymmetryModes)));
            p.addParameter('flatTopGaussianForVisualRFcenterCharacteristicRadiusEstimation', false, @islogical);
            p.addParameter('multiStartsNum', 10, @isscalar);
            p.addParameter('doDryRunFirst', false, @islogical);
            p.parse(varargin{:});


            % Generate filename for saved object
            obj.computedObjDataFileName = RetinaToVisualFieldTransformer.computedObjectDataFileName(opticsParams, targetVisualRFDoGparams);
            fprintf('Computed object will be saved to %s\n', obj.computedObjDataFileName);

            obj.theConeMosaic = theConeMosaic;
            if (isfield(theConeMosaic.coneApertureModifiers, 'shape')) && (strcmp(theConeMosaic.coneApertureModifiers.shape, 'Gaussian'))
                obj.coneCharacteristicRadiusConversionFactor = theConeMosaic.coneApertureModifiers.sigma * sqrt(2.0);
            else
                obj.coneCharacteristicRadiusConversionFactor = 0.204 * sqrt(2.0);
            end

            obj.opticsParams = opticsParams;
            obj.targetVisualRFDoGparams = targetVisualRFDoGparams;

            obj.psfWavelengthSupport = p.Results.psfWavelengthSupport;
            obj.psfCircularSymmetryMode = p.Results.psfCircularSymmetryMode;
            obj.flatTopGaussianForVisualRFcenterCharacteristicRadiusEstimation = p.Results.flatTopGaussianForVisualRFcenterCharacteristicRadiusEstimation;
            obj.multiStartsNum = p.Results.multiStartsNum;
            obj.doDryRunFirst = p.Results.doDryRunFirst;

            % Assert that the cone mosaic contains the cone pooling RF position
            % specified in the optical params struct
            coneMosaicOutline.x = theConeMosaic.eccentricityDegs(1) + theConeMosaic.sizeDegs(1)*0.5*[-1  1 1 -1 -1];
            coneMosaicOutline.y = theConeMosaic.eccentricityDegs(2) + theConeMosaic.sizeDegs(2)*0.5*[-1 -1 1  1 -1];
            in = inpolygon(obj.opticsParams.rfPositionEccDegs(1), obj.opticsParams.rfPositionEccDegs(2), ...
                coneMosaicOutline.x, coneMosaicOutline.y);
            assert(in, 'rfPosition (%2.3f,%2.3f) is outside the passed cone mosaic', ...
                obj.opticsParams.rfPositionEccDegs(1), obj.opticsParams.rfPositionEccDegs(2));

            % Assert that theConeMosaic and the optics refer to the same eye
            assert(strcmp(theConeMosaic.whichEye, obj.opticsParams.analyzedEye), ...
                'theConeMosaic and the optics params must having matching eye');

        
            % Generate the PSF to employ
            obj.vLambdaWeightedPSFandOTF();

            % Estimate the characteristic radius of the mean cone at the cone pooling RF position
            %  as projected on to visual space using the computed PSF
            dStruct = RetinaToVisualFieldTransformer.estimateConeCharacteristicRadiusInVisualSpace(...
                obj.theConeMosaic, obj.theCircularPSFData, obj.opticsParams.rfPositionEccDegs, ...
                obj.coneCharacteristicRadiusConversionFactor);
    
            if (dStruct.conesNumInRetinalPatch==0)
                error('No cones in cone mosaic')
            end

            % Compute maxSpatialSupportDegs based on desired visual RF properties
            maxSpatialSupportDegs = dStruct.visualConeCharacteristicRadiusDegs * 2.5 * ...
                   sqrt(targetVisualRFDoGparams.conesNumPooledByTheRFcenter) * ...
                   targetVisualRFDoGparams.surroundToCenterRcRatio;
            maxSpatialSupportDegs = round(100*maxSpatialSupportDegs)/100;

            
            % Crop the PSFs to maxSpatialSupportDegs to speed up computations
            obj.cropPSF(maxSpatialSupportDegs);

            % Unit weights for center cones
            indicesOfConesPooledByTheRFcenter = dStruct.indicesOfConesSortedWithDistanceToTargetRFposition(1:targetVisualRFDoGparams.conesNumPooledByTheRFcenter);
            weightsOfConesPooledByTheRFcenter = ones(1,numel(indicesOfConesPooledByTheRFcenter));

            % Compute retinal cone pooling params to generate the target visual RF
            obj.retinalRFparamsForTargetVisualRF(indicesOfConesPooledByTheRFcenter, ...
                weightsOfConesPooledByTheRFcenter, targetVisualRFDoGparams);

            % Save the computed object
            fprintf('Saving computed object to %s\n', obj.computedObjDataFileName);
            save(obj.computedObjDataFileName, 'obj');

            % Visualize the results
            obj.visualizeResults();
        end

        % Method to visualize the results
        visualizeResults(obj);

    end % Public methods

    % Private methods
    methods (Access=private)
        % Generate vLambda weighted PSF and its circularly symmetric version
        vLambdaWeightedPSFandOTF(obj);

        % Crop the PSFs
        cropPSF(obj,maxSpatialSupportDegs);

        % Method to compute the cone map for the RF center and its
        % corresponding Gaussian characteristic radius
        [visualRFcenterConeMap, visualRFcenterCharacteristicRadiusDegs, ...
         visualRFcenterCharacteristicRadiiDegs, visualRFcenterFlatTopExponents, ...
         visualRFcenterXYpos, visualRFcenterOrientationDegs] = analyzeRFcenter(obj, ...
            indicesOfConesPooledByTheRFcenter, weightsOfConesPooledByTheRFcenter, ...
            spatialSupportDegs);


        % Obtain the retinal cone pooling params (weights and indices of surround cones) by fitting the target visualRF
        dataOut = retinalRFparamsForTargetVisualRF(obj,indicesOfConesPooledByTheRFcenter, weightsOfConesPooledByTheRFcenter, targetVisualRFDoGparams);
    
    end

    % Class methods
    methods (Static)
        % Method to generate a circularly symmetric  PSF from a given PSF
        theCircularPSF = circularlySymmetricPSF(thePSF, mode);

        % Method to estimate the visually-projectected cone Rc given a target
        % position in the mosaic and corresponding PSF data
        dStruct = estimateConeCharacteristicRadiusInVisualSpace(theConeMosaic, thePSFData, theTargetPositionDegs, coneCharacteristicRadiusConversionFactor);
    
        % Method to estimate various aspects of the geometry of a 2D shape
        [theCentroid, theAxesLengths, theRotationAngle] = estimateGeometry(supportX, supportY, zData);

        % Method to analyze the effect of optics on the anatomical cone characteristic radius 
        visualConeCharacteristicRadiusDegs = analyzeVisuallyProjectedConeAperture(...
                 anatomicalConeCharacteristicRadiusDegs, thePSFData, ...
                 hFig, videoOBJ, pdfFileName);

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

        % retinal cone pooling model: arbitrary center/gaussian surround
        pooledConeIndicesAndWeights = conePoolingCoefficientsForArbitraryCenterGaussianSurround(...
            modelConstants, conePoolingParamsVector);

        % retinal cone pooling model: arbitrary center/double gaussian surround
        pooledConeIndicesAndWeights = conePoolingCoefficientsForArbitraryCenterDoubleGaussianSurround(...
            modelConstants, conePoolingParamsVector);

        % retinal cone pooling model: arbitrary center/gaussian surround
        % with arbitrary adjustments in the surround
        pooledConeIndicesAndWeights = conePoolingCoefficientsForArbitraryCenterGaussianAdjustSurround(...
            modelConstants, conePoolingParamsVector)

        % Compute the fitted visualRF from the current retinal pooling params
        [theFittedVisualRF, theRetinalRFcenterConeMap, theRetinalRFsurroundConeMap] = ...
            visualRFfromRetinalConePooling(modelConstants, retinalPoolingParams);

        % Method to compute the retinal subregion cone map by summing the
        % corresponding cone aperture maps
        retinalSubregionConeMap = retinalSubregionConeMapFromPooledConeInputs(coneRc, conePos, coneWeights, spatialSupport);

        % Method to fit a 2D Gaussian ellipsoid
        theFittedGaussian = fitGaussianEllipsoid(supportX, supportY, theRF, varargin);

        % Method to visualize the fitter param values
        visualizeFittedParamValues(retinalConePoolingParams);

        % Method to generate the datafilename where the computed object is saved
        dataFileName = computedObjectDataFileName(opticsParams, targetVisualDoGparams);
    end


end
