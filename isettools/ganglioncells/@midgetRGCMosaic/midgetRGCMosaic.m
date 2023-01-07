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

        % The positions (within the mosaic) at which we fit an RTVFT object
        theSamplingPositionGrid;

        % The # of center cones for which we fit an RTVFT object
        theConesNumPooledByTheRFcenterGrid;

        % The visual params for which we fit an RTVFT object
        theVisualSTFSurroundToCenterRcRatioGrid;
        theVisualSTFSurroundToCenterIntegratedSensitivityRatioGrid;

        % The currently employed optical image and its index in the theSamplingPositionGrid
        theCurrentOpticalImage = [];
        theCurrentOpticalImagePositionGridIndex = [];

        % Flag indicating whether this is a frozen mosaic.
        isFrozen = false;
    end % Read-only properties


    % Dependent properties
    properties (Dependent)
        % The mosaic's temporal equivalent eccentricity.
        % If the mosaic is on the temporal retinal quadrant this is: sqrt(ecc_x^2 + ecc_y^2)
        % If the mosaic is on the nasal retinal quadrant this is: sqrt((0.61 * ecc_x)^2 + ecc_y^2)
        temporalEquivalentEccentricityDegs;

        % Number of cells
        cellsNum;

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

        % Method to freeze the mosaic, removing components of the various
        %theRetinaToVisualFieldTransformerOBJs that take a lot of space
        freeze(obj);
        
        % Method to compute the response of the midgetRGCmosaic to a scene
        [midgetRGCresponses, responseTemporalSupport, ...
         noiseFreeAbsorptionsCount] = compute(obj, theScene, varargin);

        % Method to compute optics at a desired position within the mosaic
        generateOpticsAtPosition(obj, wavefrontOpticsPositionDegs);

        % Method to compute the retinal RFcenter maps - used for
        % visualization and RFoverlap analysis
        retinalRFcenterMaps = computeRetinalRFcenterMaps(obj, marginDegs, spatialSupportSamplesNum, varargin);

        % Method to generate the C/S spatial pooling RFs based on the computed list of RTVFT
        % object(s) at different positions and # of center cones
        generateCenterSurroundSpatialPoolingRFs(obj, theRetinaToVisualFieldTransformerOBJList, ...
            theSamplingPositionGrid, theConesNumPooledByTheRFcenterGrid, ...
            theVisualSTFSurroundToCenterRcRatioGrid, ...
            theVisualSTFSurroundToCenterIntegratedSensitivityRatioGrid);

        % Method to return the indices and weights of the triangulating RTVFobjects for an RGC
        [triangulatingRTVFobjIndices, triangulatingRTVFobjWeights] = ...
            triangulatingRTVFobjectIndicesAndWeights(obj,iRGC);


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

        % Method to visual the spatial RFs
        [hFig, allAxes] = visualizeSpatialRFs(obj, varargin);
        
        % Get the temporal equivalent eccentricity for either the mosaic's
        % ecc or any [Mx2] matrix of eccentricities
        eccDegs = temporalEquivalentEccentricityForEccentricity(obj, varargin);

        % Getter for dependent property cellsNum
        function val = get.cellsNum(obj)
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
        % Method to generate RetinaToVisualTransformer objects for a
        % midgetRGCMosaic (OLD)
        [RTVFTobjList, theSamplingPositionGrid, ...
         theConesNumPooledByTheRFcenterGrid, theVisualSTFSurroundToCenterRcRatioGrid, ...
         theVisualSTFSurroundToCenterIntegratedSensitivityRatioGrid] = generateRTVFobjects(...
                   ZernikeDataBase, subjectRankOrder, pupilDiameterMM, ...
                   theMidgetRGCMosaic, eccentricitySamplingGrid, ...
                   centerConnectableConeTypes, surroundConnectableConeTypes, ...
                   coneWeightsCompensateForVariationsInConeEfficiency, ...
                   CronerKaplanMultipliers, ...
                   visualRFmodel, retinalConePoolingModel, ...
                   targetSTFmatchMode, ...
                   multiStartsNum, multiStartsNumDoGFit);

        % Method to generate RetinaToVisualTransformer objects for a
        % midgetRGCMosaic (NEW)
        [RTVFTobjList, theSamplingPositionGrid, theConesNumPooledByTheRFcenterGrid, ...
          theVisualSTFSurroundToCenterRcRatioGrid, ...
          theVisualSTFSurroundToCenterIntegratedSensitivityRatioGrid] = R2VFTobjects(...
                   theMidgetRGCMosaic, ...
                   eccentricitySamplingGrid, ...
                   mosaicSurroundParams, opticsParams, fitParams);

        extraDegs = extraConeMosaicDegsForMidgetRGCSurrounds(eccentricityDegs, sizeDegs);

        % Method to generate coordinates for eccentricitySamplingGrid passed to the RTVF object
        gridCoords = eccentricitySamplingGridCoords(eccentricityDegs, sizeDegs, gridHalfSamplesNum, samplingScheme, visualizeGridCoords);

    end % Static metthods

end
