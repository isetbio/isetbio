classdef mRGCMosaic < handle
% Create an mRGCMosaic mosaic object

    % Public properties
    properties  (GetAccess=public, SetAccess=public)
        % Name of the mosaic
        name;

    end % Public properties


    % Constant properties
    properties (Constant)
        % Default optics params
        defaultOpticsParams = struct(...
            'positionDegs', [], ...
            'ZernikeDataBase', 'Polans2015', ...
            'examinedSubjectRankOrder', 6, ...
            'pupilDiameterMM', 3.0, ...
            'refractiveErrorDiopters', 0.0, ...    % use -999 for optics that do not subtract the central refraction
            'analyzedEye', 'right eye', ...
            'subjectRankingEye', 'right eye', ...
            'zeroCenterPSF', true, ...
            'psfUpsampleFactor', 1.0, ...
            'wavefrontSpatialSamples', 401); 

        defaultRetinalRFparams = struct(...
            );

        defaultVisualRFparams = struct(...
            );
        
    end % Constant properties

    % Read-only properties
    properties (GetAccess=public, SetAccess=private)

        % The source lattice specification
        sourceLatticeSizeDegs;

        % Eccentricity [x,y] of the center of the mRGC mosaic (degs)
        eccentricityDegs;

        % Eccentricity [x,y] of the center of the mRGC mosaic (microns)
        eccentricityMicrons;

        % Size of the mRGC mosaic [width, height]
        sizeDegs;

        % Either 'right eye' or 'left eye'
        whichEye;

        % The input cone mosaic (@cMosaic)
        inputConeMosaic;

        % Extra size of the input cone mosaic (to allow for cone inputs to
        % the RF surrounds)
        extraDegsForInputConeMosaic;

        % [0: minimize chromatic variance 1: minimize spatial variance]
        chromaticSpatialVarianceTradeoff;

        % [n x 2] matrix of rgc rf positions, in microns & degrees
        rgcRFpositionsMicrons;
        rgcRFpositionsDegs;

        % [n x 1] vector of rgc rf spacings, in microns & degrees
        rgcRFspacingsMicrons;
        rgcRFspacingsDegs;

        % The factor by which we convert nasal eccentricities to their
        % temporal equivalents
        temporalEquivantEccentricityFactor;

        % Computed params
        % Sparse [conesNum x rgcRFsNum]  connectivity matrix 
        % To find which cones are connected to an rgcRF:
        %  connectivityVector = full(squeeze(rgcRFcenterConeConnectivityMatrix(:, rgcRF)));
        %  inputConeIndices = find(abs(connectivityVector) > 0.0001);
        rgcRFcenterConeConnectivityMatrix = [];

        rgcRFcenterConePoolingMatrix = [];
        rgcRFsurroundConePoolingMatrix = [];

        % The optics used to optimize surround cone weights
        % so as to generate RFs with the target visual properties
        theNativeOptics;

        % The params used to generate theNativeOptics
        theNativeOpticsParams;

        % Encodes what generation stage the mosaic is in.
        generationStage;

    end % Read-only properties

    % Dependent properties
    properties (Dependent)
        % The mosaic's temporal equivalent eccentricity.
        % If the mosaic is on the temporal retinal quadrant this is: sqrt(ecc_x^2 + ecc_y^2)
        % If the mosaic is on the nasal retinal quadrant this is: sqrt((0.61 * ecc_x)^2 + ecc_y^2)
        temporalEquivalentEccentricityDegs;

        % Number of cells
        rgcsNum;

        % Retinal quadrant (i.e. nasal, temporal, inferior, superior)
        horizontalRetinalMeridian;
        verticalRetinalMeridian;
    end

    % Private properties
    properties (GetAccess=private, SetAccess=private)
        % The MosaicConnectorOBJ used to connect cones to midget RGCs
        theMosaicConnectorOBJ;

        % Cache for various visualizations
        visualizationCache = [];
    end

    % Public methods
    methods
        
        % Constructor
        function obj = mRGCMosaic(varargin)

            % Parse input
            p = inputParser;
            p.addParameter('name', 'mRGC mosaic', @ischar);
            p.addParameter('sourceLatticeSizeDegs', 58, @isscalar);
            p.addParameter('whichEye', 'right eye', @(x)(ischar(x) && (ismember(x, {'left eye', 'right eye'}))));
            
            % Eccentricity/Position: if positionDegs is sent in, it overwrites eccentricityDegs
            p.addParameter('eccentricityDegs', [], @(x)(isnumeric(x) && (numel(x) == 2)));
            p.addParameter('positionDegs', [], @(x)(isnumeric(x) && (numel(x) == 2)));
            p.addParameter('sizeDegs', [0.4 0.4], @(x)(isnumeric(x) && (numel(x) == 2)));
            p.addParameter('temporalEquivantEccentricityFactor', 'ISETBioMosaicsBased', @(x)(ischar(x) && (ismember(x, {'ISETBioMosaicsBased', 'WatanabeRodieckBased'}))));
            
            % Custom Degs <-> MM conversion functions
            p.addParameter('customDegsToMMsConversionFunction', @(x)RGCmodels.Watson.convert.rhoDegsToMMs(x), @(x) (isempty(x) || isa(x,'function_handle')));
            p.addParameter('customMMsToDegsConversionFunction', @(x)RGCmodels.Watson.convert.rhoMMsToDegs(x), @(x) (isempty(x) || isa(x,'function_handle')));
            
            p.addParameter('chromaticSpatialVarianceTradeoff', 1.0, @(x)(isscalar(x)&&((x>=0)&&(x<=1))));
            p.parse(varargin{:});

            obj.name = p.Results.name;
            obj.sourceLatticeSizeDegs = p.Results.sourceLatticeSizeDegs;
            obj.chromaticSpatialVarianceTradeoff = p.Results.chromaticSpatialVarianceTradeoff;
            obj.temporalEquivantEccentricityFactor = p.Results.temporalEquivantEccentricityFactor;
            
            % Generate the input cone mosaic. This also sets various
            % other properties that are matched to the input cone
            % mosaic, such as:
            % - eccentricityDegs
            % - sizeDegs
            % - whichEye
            obj.generateInputConeMosaic(p.Results);

            % Set the eccentricity in microns
            obj.eccentricityMicrons = obj.inputConeMosaic.eccentricityMicrons;

            % Wire the RF centers. This configures and runs theobj.theMosaicConnectorOBJ
            % which wires cones to midget RGC RF center (no overlap)
            obj.generateRFpositionsAndWireTheirCenters();

            % Eliminate RGCs near the border
            obj.cropRGCsOnTheBorder();
        end % Constructor

        % Method to generate the native optics for the mosaic.
        % These optics form the basis on which surround cone weights are optimized
        % for wiring so as to generate desired visual RF properties
        generateNativeOptics(obj, opticsParams);

        % Method to visualize the mRGCmosaic and its activation
        visualize(obj, varargin);

        % Method to bake in center/surround cone pooling weights.
        % This method replaces the rgcRFcenterConeConnectivityMatrix with
        % the rgcRFcenterConePoolingMatrix & rgcRFsurroundConePoolingMatrix
        % and can only be called by the MosaicPoolingOptimizer which
        % optimized these weights based on input cone mosaic STF responses
        % computed by computeInputConeMosaicSTFresponsesForNativeOptics
        bakeInConePoolingMatrices(obj, centerConePoolingMatrix, surroundConePoolingMatrix);

        % Getter for dependent property cellsNum
        function val = get.rgcsNum(obj)
            val = size(obj.rgcRFpositionsDegs,1);
        end

        % Getter for dependent property temporalEquivalentEccentricityDegs
        function val = get.temporalEquivalentEccentricityDegs(obj)
            val = obj.temporalEquivalentEccentricityForEccentricity(obj.eccentricityDegs);
        end

        % Getter for dependent property horizontalRetinalMeridian
        function val = get.horizontalRetinalMeridian(obj)
            val = obj.inputConeMosaic.horizontalRetinalMeridian;
        end

        % Getter for dependent property verticalRetinalMeridian
        function val = get.verticalRetinalMeridian(obj)
            val = obj.inputConeMosaic.verticalRetinalMeridian;
        end

    end % Public methods

    % Private methods
    methods (Access=private)
        generateInputConeMosaic(obj, pResults);
        generateRFpositionsAndWireTheirCenters(obj);
        
        % Method to crop RGCs on the border
        cropRGCsOnTheBorder(obj);
    end % Private methods


    % Static methods
    methods (Static)
    end % Static methods

end



 