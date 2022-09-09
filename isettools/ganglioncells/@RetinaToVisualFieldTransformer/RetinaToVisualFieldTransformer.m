classdef RetinaToVisualFieldTransformer < handle

    % Public properties (read-only)
    properties (SetAccess = private)

       % The wavelength support for the employed PSF
       psfWavelengthSupport;

       % The PSF circular symmetry mode
       psfCircularSymmetryMode;

       % The @cMosaic object
       theConeMosaic;

       % Various optics params
       opticsParams;

       % Subject ID for the chosen database and ranking
       testSubjectID;

       % Whether the subject requires subtraction of central refraction
       subtractCentralRefraction;

       % The computed (vLambda weighted, or not) PSFs
       thePSFData;
       theCircularPSFData;
    end

    % Constant properties
    properties (Constant, Hidden)
        Artal = 'Artal2012';
        Polans = 'Polans2015';
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
        function obj = RetinaToVisualFieldTransformer(theConeMosaic, opticsParams, targetVisualRFDoGparams, varargin) 
            % Parse input
            p = inputParser;
            p.addParameter('psfWavelengthSupport', [], @(x)(isvector(x)));
            p.addParameter('psfCircularSymmetryMode', ...
                RetinaToVisualFieldTransformer.psfCircularSymmetryModeBestResolution, ...
                @(x)(ismember(x, RetinaToVisualFieldTransformer.validPSFcircularSymmetryModes)));
            p.parse(varargin{:});

            obj.theConeMosaic = theConeMosaic;
            obj.opticsParams = opticsParams;
            
            obj.psfWavelengthSupport = p.Results.psfWavelengthSupport;
            obj.psfCircularSymmetryMode = p.Results.psfCircularSymmetryMode;

            % Assert that the cone mosaic contains the cone pooling RF position
            % specified in the optical params struct
            coneMosaicOutline.x = theConeMosaic.eccentricityDegs(1) + theConeMosaic.sizeDegs(1)*0.5*[-1  1 1 -1 -1];
            coneMosaicOutline.y = theConeMosaic.eccentricityDegs(2) + theConeMosaic.sizeDegs(2)*0.5*[-1 -1 1  1 -1];
            [in,on] = inpolygon(obj.opticsParams.rfPositionEccDegs(1), obj.opticsParams.rfPositionEccDegs(2), ...
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
                obj.theConeMosaic, obj.theCircularPSFData, obj.opticsParams.rfPositionEccDegs);
    
            if (dStruct.conesNumInRetinalPatch==0)
                error('No cones in cone mosaic')
            end

            % Compute maxSpatialSupportDegs based on desired visual RF properties
            maxSpatialSupportDegs = mean([dStruct.visualConeCharacteristicRadiusDegs dStruct.anatomicalConeCharacteristicRadiusDegs]) * 2.0 * ...
                   sqrt(targetVisualRFDoGparams.conesNumPooledByTheRFcenter) * ...
                   targetVisualRFDoGparams.surroundToCenterRcRatio;
            maxSpatialSupportDegs = round(100*maxSpatialSupportDegs)/100;

            
            % Crop the PSFs to maxSpatialSupportDegs to speed up computations
            obj.cropPSF(maxSpatialSupportDegs);

            % Compute retinal cone pooling params to generate the target
            % visual RFs

            STOPED HERE
            
           [retinalRFparams, weightsComputeFunctionHandle, ...
            targetVisualRF, spatialSupportDegs, modelConstants] = obj.retinalRFparamsForTargetVisualRF();

        end

    end % Public methods

    % Private methods
    methods (Access=private)
        % Generate vLambda weighted PSF and its circularly symmetric version
        vLambdaWeightedPSFandOTF(obj);

        % Crop the PSFs
        cropPSF(obj,maxSpatialSupportDegs);
    end

    % Class methods
    methods (Static)
        % Method to generate a circularly symmetric  PSF from a given PSF
        theCircularPSF = circularlySymmetricPSF(thePSF, mode);

        % Method to estimate the visually-projectected cone Rc given a target
        % position in the mosaic and corresponding PSF data
        dStruct = estimateConeCharacteristicRadiusInVisualSpace(theConeMosaic, thePSFData, theTargetPositionDegs);
    
        % Method to estimate various aspects of the geometry of a 2D shape
        [theCentroid, theAxesLengths, theRotationAngle] = estimateGeometry(supportX, supportY, zData);

        % Method to analyze the effect of optics on the anatomical cone characteristic radius 
        visualConeCharacteristicRadiusDegs = analyzeVisuallyProjectedConeAperture(...
                 anatomicalConeCharacteristicRadiusDegs, thePSFData, ...
                 hFig, videoOBJ, pdfFileName);
        
        % Method to fit a 2D Gaussian ellipsoid
        theFittedGaussian = fitGaussianEllipsoid(supportX, supportY, theRF, varargin);
    end


end
