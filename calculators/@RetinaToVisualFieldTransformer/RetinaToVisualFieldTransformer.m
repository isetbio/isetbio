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
        temporalRetinaQuadrant = 'temporal retina';
        nasalRetinaQuadrant = 'nasal retina';
        inferiorRetinaQuadrant = 'inferior retina';
        superiorRetinaQuadrant = 'superior retina';
        leftEye = 'left eye';
        rightEye = 'right eye';
    end

    properties (Constant)
        
        validZernikeDataBases = {...
            RetinaToVisualFieldTransformer.Artal ...
            RetinaToVisualFieldTransformer.Polans ...
            };

        validRetinalQuadrants = { ...
            RetinaToVisualFieldTransformer.temporalRetinaQuadrant ...
            RetinaToVisualFieldTransformer.nasalRetinaQuadrant ...
            RetinaToVisualFieldTransformer.inferiorRetinaQuadrant ...
            RetinaToVisualFieldTransformer.superiorRetinaQuadrant ...
            };

        validEyes = {...
            RetinaToVisualFieldTransformer.leftEye ...
            RetinaToVisualFieldTransformer.rightEye ...
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

        % Estimate the characteristic radius of the mean cone at the center of a
        % cone mosaic centered at a given (x,y) eccentricity in degrees
        % as projected on to visual space using optics at the same eccentricity
        dStruct = estimateConeCharacteristicRadiusInVisualSpace(obj,...
            whichEye, eccDegs, testSubjectID, pupilDiameterMM, ...
            dataFileName, varargin);


    end % Public methods

    % Private methods
    methods (Access=private)
    end

    % Class methods
    methods (Static)
        % Method to return the range of horizontal and vertical
        % eccentricites for a retinal quadrant and eye
        [horizontalEccDegs, verticalEccDegs, eccDegsForPlotting] = ...
                eccentricitiesForQuadrant(retinalQuadrant, whichEye, maxEcc);

    end


end
