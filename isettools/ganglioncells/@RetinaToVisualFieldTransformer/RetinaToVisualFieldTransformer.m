classdef RetinaToVisualFieldTransformer

    % Public properties (read-only)
    properties (SetAccess = private)
       % The Zernike coefficients to use
       ZernikeDataBase;
    end

    % Constant properties related to figure generation
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
        function obj = RetinaToVisualFieldTransformer(varargin) 
            % Parse input
            p = inputParser;
            p.addParameter('ZernikeDataBase', 'Artal2012', @(x)(ismember(x, RetinaToVisualFieldTransformer.validZernikeDataBases)));
            p.parse(varargin{:});
            obj.ZernikeDataBase = p.Results.ZernikeDataBase;
        end

        % Return the ID of a subject with a specific rank order for a given eye
        testSubjectID = subjectWithRankInEye(obj, subjectRankOrder, retinalQuadrant, whichEye);

        % Generate vLambda-weighted PSF
        [thePSFData, theCircularPSFData] = vLambdaWeightedPSFandOTF(obj, cm, testSubjectID, ...
            pupilDiameterMM, wavefrontSpatialSamples, maxSpatialSupportDegs, circularSymmetryGenerationMode);

        % Estimate the characteristic radius of the mean cone at the center of a
        % cone mosaic centered at a given (x,y) eccentricity in degrees
        % as projected on to visual space using optics at the same eccentricity
        dStruct = estimateConeCharacteristicRadiusInVisualSpace(obj,...
            whichEye, eccDegs, testSubjectID, pupilDiameterMM, ...
            dataFileName, varargin);

        % Deconvole the visualRF based on V-lambda weighted optics at a
        % given eccentricity, subject, and pupil size to extrace the
        % retinal cone pooling weights for the center/surround mechanisms
        [retinalRFparams, weightsComputeFunctionHandle, ...
         targetVisualRF, spatialSupportDegs, theCircularPSFData] = retinalRFparamsForTargetVisualRF(obj, ...
          visualRFDoGparams, eccDegs, ...
          subjectEye, subjectID, subjectPupilDiameterMM, varargin);

    end % Public methods

    % Private methods
    methods (Access=private)
    end

    % Class methods
    methods (Static)
        % Method to return the centroid (in pixels) and in 
        % units of the support and the rotation in degrees of the zData
        [theCentroid, RcX, RcY, theRotationAngle] = ...
        estimateGeometry(supportX, supportY, zData);
        
        theCircularPSF = circularlySymmetricPSF(thePSF, mode);

        data = centerAndRotatePSF(data);

        visualConeCharacteristicRadiusDegs = analyzeEffectOfPSFonConeAperture(...
            anatomicalConeCharacteristicRadiusDegs, thePSFData, ...
            hFig, videoOBJ, pdfFileName);

        [RcDegs, rotationDegs, flatTopGaussianExponent, ...
         visualRFcenterConeMap, retinalRFcenterConeMap, ...
         anatomicalConeCharacteristicRadiusDegs] = estimateVisualRcFromNumberOfConesInRFcenter(cm, ...
            conesNumPooledByTheRFcenter, theCircularPSFData, spatialSupportDegs);

        [RFcenterRcDegs, RF2D] = computeRetinalRFRcDegsFromItsPooledConeInputs(coneRcDegs, ...
            conePosDegs, spatialSupportDegs);
    
        [theFittedGaussianCharacteristicRadiusDegs, visualConeCharacteristicMinorMajorRadiiDegs, ...
         theFittedGaussianRotationDegs, theFittedGaussianFlatTopExponent, ...
         theFitted2DGaussian, XYcenter, XRange, YRange] = ...
            fitGaussianEllipsoid(supportX, supportY, theRF, varargin);


        [theFittedGaussianCharacteristicRadiusDegs, theFittedGaussianEllpsoid, ...
            XYcenter, XRange, YRange, theFittedGaussianCharacteristicRadiiDegs] = ...
            fitGaussianToPooledConeApertures(supportX, supportY, theRF, initialParams, lowerBounds, upperBounds);
        RF2F = gaussian2D(params,xydata);

        paramsVector = paramsStructToParamsVector(paramsStruct, retinalConePoolingModel);
        paramsStruct = paramsVectorToParamsStruct(paramsVector, retinalConePoolingModel);

        [theFittedRFparamsStruct, theFittedRF] = fitDoGModelToRF(spatialSupportDegs, RF, rotationDegs, RcDegs);
        RF2D = diffOfGaussiansRF(p, spatialSupportDegs);

        [theFittedRFparamsStruct, theFittedRF] = fitDiffOfGaussianCenterAndDoubleExponentSurroundModelToRF(spatialSupportDegs, RF, rotationDegs, RcDegs);
        RF2D = diffOfGaussianCenterAndDoubleExponentSurround(p, spatialSupportDegs);

        [retinalRFparams, weightsComputeFunctionHandle, theFittedVisualRF, ...
         theFittedRetinalRFcenter, theFittedRetinalRFsurround] = fitVisualRFByAdjustingRetinalPoolingParameters(...
                modelConstants, theVisualRF)

        pooledConeIndicesAndWeightsStruct = retinalConeWeightsFromDoGDEmodelParameters(retinalRFDoGDEparams, ...
            conesNumPooledByTheRFcenter, spatialSupportX, spatialSupportY, cm);

        pooledConeIndicesAndWeightsStruct = retinalConeWeightsFromDoGmodelParametersForTargetRFCenterCones(...
           retinalRFDoGparams, cm, targetRFCenterConesIndices);

        pooledConeIndicesAndWeightsStruct = retinalConeWeightsFromDoGmodelParameters(retinalRFDoGparams, ...
            conesNumPooledByTheRFcenter, spatialSupportX, spatialSupportY, cm);

        retinalRF = computeRetinalRFfromWeightsAndApertures(consideredConeIndices, coneIndices, coneWeights, ...
            coneApertureDiametersDegs, coneRFpositionsDegs, Xdegs, Ydegs, RFcenter, minConeWeight);

        [retinalRFcenter2D, retinalRFsurround2D] = generateRFsubregionMapsFromPooledCones(...
            spatialSupportX, spatialSupportY, cm, pooledConeIndicesAndWeightsStruct);

        generateFigure(thePSFData, theCircularPSFData, RF2DData, ...
            visualRFparams, retinalRFparams, achievedVisualRFparams, ...
            eccDegs, testSubjectID, maxSpatialSupportDegs, rfSupportX, rfSupportY, figNo);

    end


end
