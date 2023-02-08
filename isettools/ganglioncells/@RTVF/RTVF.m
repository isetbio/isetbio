classdef RTVF < handle


    % Public properties (read-only)
    properties (SetAccess = private)
       % The @cMosaic object
       theConeMosaic;

       % The input optics params
       opticsParams;

       % The target visual params
       targetVisualRFDoGparams;

       % The exports directory
       exportsDirectory;

       % L- and M-cone weighted PSFs 
       % appropriate for the optical
       % position 
       theSpectrallyWeightedPSFData;

       % Visualization optiocs
       visualizeSpectrallyWeightedPSFs;

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
            'visualRFmodel', 'arbitraryShapedCenter_GaussianShapedSurround', ... 
            'retinalConePoolingModel', 'arbitraryCenterConeWeights_doubleExpH1cellIndex1SurroundWeights', ...
            'centerConnectableConeTypes', [cMosaic.LCONE_ID cMosaic.MCONE_ID], ...
            'surroundConnectableConeTypes', [cMosaic.LCONE_ID cMosaic.MCONE_ID], ...
            'coneWeightsCompensateForVariationsInConeEfficiency', true,  ...
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
            p.addParameter('exportsDirectory', '', @(x)(isempty(x)||ischar(x)));
            p.addParameter('visualizeSpectrallyWeightedPSFs', false, @islogical);
            p.parse(varargin{:});

            % Handle optional inputs
            obj.exportsDirectory = p.Results.exportsDirectory;
            obj.visualizeSpectrallyWeightedPSFs = p.Results.visualizeSpectrallyWeightedPSFs;
       
            % Handle empty cone mosaic
            if (isempty(theConeMosaic))
                % One L- and one M-cone only
                obj.theConeMosaic = cMosaic('sizeDegs', [0.02 0.01]);
                obj.theConeMosaic.reassignTypeOfCones([], cMosaic.KCONE_ID);
                obj.theConeMosaic.reassignTypeOfCones(1, cMosaic.LCONE_ID);
                obj.theConeMosaic.reassignTypeOfCones(2, cMosaic.MCONE_ID);
                obj.theConeMosaic.visualize()
            else
                obj.theConeMosaic = theConeMosaic;
            end

            % Handle empty opticsParams
            if (isempty(theOpticsParams))
                obj.opticsParams = RTVF.defaultOpticsParams;
            else
                obj.opticsParams = theOpticsParams;
            end

            % Handle empty argetVisualRFDoGparams
            if (isempty(theTargetVisualRFDoGparams))
                obj.targetVisualRFDoGparams = RTVF.defaultTargetVisualRFDoGParams;
            else
                obj.targetVisualRFDoGparams = theTargetVisualRFDoGparams;
            end
            

            % Compute spectrally weighted PSFs.
            % This sets the obj.theSpectrallyWeightedPSFData parameter
            % and updates the obj.opticsParams
            obj.spectrallyWeightedPSFs();

        end % 
        % Constructor
    end % public methods

    % Class methods
    methods (Static)
    end

end

