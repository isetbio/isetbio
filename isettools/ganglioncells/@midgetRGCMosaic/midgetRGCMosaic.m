classdef midgetRGCMosaic < handle
% Create a midgetRGCMosaic mosaic object

    % Public properties
    properties  (GetAccess=public, SetAccess=public)
        % Name of the mosaic
        name;

    end % public properties

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

        % Overlap of RFs (equivalent to 1/normalizedNearestNeighorDistance)
        % See: Gauthier et al., 2009, "Uniform Signal Redundancy in Primate Retina" 
        rfOverlapRatio;

        % [0: minimize chromatic variance 1: minimize spatial variance]
        chromaticSpatialVarianceTradeoff;

        % Computed params
        % Sparse [conesNum x rgcRFsNum]  connectivity matrix 
        % To find which cones are connected to a n rgcRF:
        %  connectivityVector = full(squeeze(rgcRFcenterConeConnectivityMatrix(:, rgcRF)));
        %  inputConeIndices = find(abs(connectivityVector) > 0.0001);
        rgcRFcenterConeConnectivityMatrix;

        % The conePoolingMatrices are sparse [conesNum x rgcRFsNum] matrices 
        % computed based on the information in theRetinaToVisualFieldTransformerOBJList
        rgcRFcenterConePoolingMatrix;
        rgcRFsurroundConePoolingMatrix;

        % [n x 2] matrix of rgc rf positions, in microns & degrees
        rgcRFpositionsMicrons;
        rgcRFpositionsDegs;

        % [n x 1] vector of rgc rf spacings, in microns & degrees
        rgcRFspacingsMicrons;
        rgcRFspacingsDegs;

        % The factor by which we convert nasal eccentricities to their
        % temporal equivalents
        temporalEquivantEccentricityFactor;

        % The MosaicConnectorOBJ used to connect cones to midget RGCs
        theMosaicConnectorOBJ;

        % The RetinalToVisualFieldTransformerList cell array which contains 
        % RTVFT objects for a grid of positions, conesNumInRFcenter,
        % and visual STF properties:
        % - optics (position) 
        % - conesNumInRFcenter
        % - visual STF properties (C/S radius, C/S int.ratio)
        theRetinaToVisualFieldTransformerOBJList;

        theOpticsPositionGrid;
        theConesNumPooledByTheRFcenterGrid;
        theVisualSTFSurroundToCenterRcRatioGrid;
        theVisualSTFSurroundToCenterIntegratedSensitivityRatioGrid;

    end % Read-only properties


    % Dependent properties
    properties (Dependent)
        % The mosaic's temporal equivalent eccentricity.
        % If the mosaic is on the temporal retinal quadrant this is: sqrt(ecc_x^2 + ecc_y^2)
        % If the mosaic is on the nasal retinal quadrant this is: sqrt((0.61 * ecc_x)^2 + ecc_y^2)
        temporalEquivalentEccentricityDegs;

        % Retinal quadrant (i.e. nasal, temporal, inferior, superior)
        horizontalRetinalMeridian;
        verticalRetinalMeridian;
    end


    % Private properties
    properties (GetAccess=private, SetAccess=private)
    
    end


    % Public methods
    methods
        
        % Constructor
        function obj = midgetRGCMosaic(varargin)
            % Parse input
            p = inputParser;
            p.addParameter('name', 'mRGC mosaic', @ischar);
            p.addParameter('inputConeMosaic', [], @(x)(isempty(x) || isa(x, 'cMosaic')));
           
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
            
            % RF overlap: 0.5 results in a median NNND of around 2.0
            p.addParameter('rfOverlapRatio', 0.0, @(x)(isscalar(x)&&((x>=0)&&(x<=1))));
            p.addParameter('chromaticSpatialVarianceTradeoff', 1.0, @(x)(isscalar(x)&&((x>=0)&&(x<=1))));
            p.parse(varargin{:});

            obj.name = p.Results.name;
            obj.sourceLatticeSizeDegs = p.Results.sourceLatticeSizeDegs;
            obj.chromaticSpatialVarianceTradeoff = p.Results.chromaticSpatialVarianceTradeoff;
            obj.temporalEquivantEccentricityFactor = p.Results.temporalEquivantEccentricityFactor;

            % Deal with the input cone mosaic. If an input cMosaic is not specified we generate one here
            if (~isempty(p.Results.inputConeMosaic))
                fprintf(2,'An input cone mosaic was passed: ignoring set key/value pairs for properties inherited by the cMosaic.\n')
                pause
            else
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

                % Adjust the RF overlap. This sets the following properties:
                % - rgcRFcenterConeConnectivityMatrix
                % - rgcRFpositionsMicrons, rgcRFpositionsDegs
                % - rgcRFspacingsMicrons, rgcRFspacingsDegs
                obj.adjustRFoverlap(p.Results.rfOverlapRatio);
            end
        end % Constructor

        % Method to compute the response of the midgetRGCmosaic to a scene
        [responses, responseTemporalSupport] = compute(obj, theScene, varargin);

        % Method to compute the retinal RFcenter maps - used for
        % visualization and RFoverlap analysis
        retinalRFcenterMaps = computeRetinalRFcenterMaps(obj, marginDegs, spatialSupportSamplesNum, varargin);

        % Method to generate the C/S spatial pooling RF based on a passed RTVFT
        % object(s), which encodes C/S weights to achieve a desired visual STF
        % for specific optics and # of center cones
        generateCenterSurroundSpatialPoolingRF(obj, theRetinaToVisualFieldTransformerOBJList, ...
            theOpticsPositionGrid, theConesNumPooledByTheRFcenterGrid, ...
            theVisualSTFSurroundToCenterRcRatioGrid, ...
            theVisualSTFSurroundToCenterIntegratedSensitivityRatioGrid);

        % Method to adjust the midgetRGC RF overlap
        adjustRFoverlap(obj, overlapRatio);

        % Method to analyze the stats of the achieved RF overlap and return
        % the NormalizedNearestNeighborDistances between neighboring pairs of RFs
        % as per Gauthier, Chichilinsky et al (2009)
        [NNNDs, NNNDtuplets, RGCdistances, distancesFromMosaicCenterDegs, targetRGCindices] = analyzeRetinalRFoverlap(obj, varargin);

        % Method to visualize the mosaic
        [figureHandle, axesHandle] = visualize(obj, varargin);
        
        % Method to visualize the retinal RF of a single RGC
        visualizeSingleRetinalRF(obj,theRGCindex, varargin);

        % Method to visualize aspects of the retinal RFs
        visualizeRetinalRFs(obj, varargin);

        % Visualize the connectivity of the RF center to the cones
        visualizeRFcenterConnectivity(obj, varargin);

        % Get the temporal equivalent eccentricity for either the mosaic's
        % ecc or any [Mx2] matrix of eccentricities
        eccDegs = temporalEquivalentEccentricityForEccentricity(obj, varargin);

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
        extraDegs = extraConeMosaicDegsForMidgetRGCSurrounds(eccentricityDegs, sizeDegs);
    end % Static metthods

end
