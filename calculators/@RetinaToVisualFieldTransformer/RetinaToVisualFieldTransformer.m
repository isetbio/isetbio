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

        [theOTFData, thePSFData] = vLambdaWeightedPSFandOTF(obj, cm, testSubjectID, pupilDiameterMM);

        % Estimate the characteristic radius of the mean cone at the center of a
        % cone mosaic centered at a given (x,y) eccentricity in degrees
        % as projected on to visual space using optics at the same eccentricity
        dStruct = estimateConeCharacteristicRadiusInVisualSpace(obj,...
            whichEye, eccDegs, testSubjectID, pupilDiameterMM, ...
            dataFileName, varargin);

        % Deconvole the visualRF based on V-lambda weighted optics at a
        % given eccentricity, subject, and pupil size
        deconvolvedVisualRFdata = deconvolveVisualRFgivenOptics(obj, ...
            visualRFdata, conesNumPooledByTheRFcenter, eccDegs, ...
            whichEye, testSubjectID, pupilDiameterMM, ...
            regularizationAlpha, varargin);

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

        [theFittedGaussianCharacteristicRadiusDegs, visualConeCharacteristicMinorMajorRadiiDegs, ...
            theFitted2DGaussian, XYcenter, XRange, YRange] = ...
            fitGaussianEllipsoid(supportX, supportY, theRF);

    end


end
