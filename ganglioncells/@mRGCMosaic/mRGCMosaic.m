classdef mRGCMosaic < handle
% Create an mRGCMosaic mosaic object

    % Public properties
    properties  (GetAccess=public, SetAccess=public)
        % Name of the mosaic
        name;

        % Noise flag for mRGC responses
        noiseFlag;

        % Random seed
        randomSeed = [];

        % vMembrane noise sigma
        vMembraneGaussianNoiseSigma = 0.07;
    end % Public properties


    % Constant properties
    properties (Constant)

        % Valid noise flags
        validNoiseFlags = {'none', 'frozen', 'random'};

        % Default optics params
        defaultOpticsParams = struct(...
            'positionDegs', [], ...
            'ZernikeDataBase', 'Polans2015', ...
            'examinedSubjectRankOrder', 6, ...
            'pupilDiameterMM', 3.0, ...
            'refractiveErrorDiopters', 0.0, ...    % use -999 for optics that do not subtract the central refraction
            'noLCA', false, ...
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

        % Sparse [conesNum x rgcRFsNum] cone pooling matrices for RF center
        % and RF surrounds.
        % To find which cones are connected to the surround of an rgcRF:
        %  connectivityVector = full(squeeze(rgcRFsurroundConePoolingMatrix(:, rgcRF)));
        %  surroundConeIndices = find(abs(connectivityVector) > 0.0001);
        rgcRFcenterConePoolingMatrix = [];
        rgcRFsurroundConePoolingMatrix = [];

        % Peak gain for each RGC. 
        rgcRFgains = [];

        % The optics used to optimize surround cone weights
        % so as to generate RFs with the target visual properties
        % These are the optics at the mosaic'c center
        theNativeOptics;

        % The params used to generate theNativeOptics
        theNativeOpticsParams;

        % Custom optics. These are optics generated at a position other
        % that the center of the mosaic
        theCustomOptics;

        % The params used to generate theCustomOptics
        theCustomOpticsParams;

        % Encodes what generation stage the mosaic is in.
        generationStage;

        % Cell array of meta data structs at each intermediate connectivity stage of theMosaicConnectorOBJ 
        intermediateMetaDataStructs;

        % Cell array of figure handles for each intermediate connectivity stage of theMosaicConnectorOBJ 
        mosaicConnectorIntermediateFigureHandles;

        % Cache for various visualizations
        visualizationCache = [];
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
            p.addParameter('withInputConeMosaic', [],  @(x) (isempty(x) || isa(x, 'cMosaic')));
            p.addParameter('temporalEquivantEccentricityFactor', 'ISETBioMosaicsBased', @(x)(ischar(x) && (ismember(x, {'ISETBioMosaicsBased', 'WatanabeRodieckBased'}))));
            
            % Custom Degs <-> MM conversion functions
            p.addParameter('customDegsToMMsConversionFunction', @(x)RGCmodels.Watson.convert.rhoDegsToMMs(x), @(x) (isempty(x) || isa(x,'function_handle')));
            p.addParameter('customMMsToDegsConversionFunction', @(x)RGCmodels.Watson.convert.rhoMMsToDegs(x), @(x) (isempty(x) || isa(x,'function_handle')));
            
            p.addParameter('maxConeInputsPerRGCToConsiderTransferToNearbyRGCs', MosaicConnector.maxSourceInputsToConsiderTransferToNearbyDestinationRF, @isscalar);
            p.addParameter('maxMeanConeInputsPerRGCToConsiderSwappingWithNearbyRGCs', MosaicConnector.maxMeanSourceInputsToConsiderSwappingWithNearbyDestinationRF, @isscalar);
            p.addParameter('chromaticSpatialVarianceTradeoff', 1.0, @(x)(isscalar(x)&&((x>=0)&&(x<=1))));
            p.addParameter('visualizeIntermediateStagesOfCenterConnectivityOptimization', false, @islogical);
            p.addParameter('saveIntermediateStagesOfCenterConnectivityOptimization', false, @islogical);
            
            % Compute-ready mosaic params
            p.addParameter('noiseFlag', 'random');
            p.addParameter('randomSeed', [], @(x)(isempty(x) || isscalar(x)));
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
            obj.generateRFpositionsAndWireTheirCenters(p.Results);

            % Eliminate RGCs near the border
            obj.cropRGCsOnTheBorder();

            % Compute-ready params
            obj.noiseFlag = p.Results.noiseFlag;
        end % Constructor

        % Method to crop a compute-ready mRGCMosaic to size at an eccentricity
        cropToSizeAtEccentricity(obj, sizeDegs, eccentricityDegs, varargin);

        % Method to set the peak gain of each mRGC 
        setPeakGains(obj, method, methodParams);

        % Compute method
        [noiseFreeMRGCresponses, noisyMRGCresponseInstances, responseTemporalSupportSeconds] = compute(obj, ...
            theConeMosaicResponse, theConeMosaicResponseTemporalSupportSeconds, varargin);

        % Method to compute noisy mRGCresponse instances
        noisyMRGCresponseInstances = noisyInstances(obj, noiseFreeMRGCresponses, varargin);

        % Method to generate the visualization cache
        generateVisualizationCache(obj, xSupport, ySupport, centerSubregionContourSamples, contourGenerationMethod);

        % Method to visualize the mRGCmosaic and its activation
        visualize(obj, varargin);

        % Method to visualize the spatial RF of an RGC near a target
        % positon, with a target # of center cones, and a target center
        % cone majority type
        theVisualizedRGCindex = visualizeRetinalConePoolingRFmapNearPosition(obj, ...
            targetRGCposition, targetCenterConesNum, ...
            targetCenterConeMajorityType, varargin);

        % Method to visualize the spatial RF of an RGC with a specific index
        visualizeRetinalConePoolingRFmapOfRGCwithIndex(obj,theRGCindex, varargin);

        % Method to generate optics for the mosaic.
        dataOut = generateOptics(obj, opticsParams);
        
        % Method to set either the native (apropriate for the mosaic's center)
        % or the custom optics (at an arbitrary position within the mosaic)
        setTheOptics(obj, opticsParams);

        % Method to visualize optics at a number of eccentricities
        hFig = visualizeOpticsAtEccentricities(obj, eccDegs, opticsParams, tickSeparationArcMin);

        % Method to bake in center/surround cone pooling weights.
        % This method replaces the rgcRFcenterConeConnectivityMatrix with
        % the rgcRFcenterConePoolingMatrix & rgcRFsurroundConePoolingMatrix
        % It is called by the MosaicPoolingOptimizer.generateComputeReadyMidgetRGCMosaic()
        % method which computes optimal center/surround weights for all RGCs in the mosaic
        % by fitting the mRGC model to the Croner&Kaplan STF data 
        bakeInConePoolingMatrices(obj, centerConePoolingMatrix, surroundConePoolingMatrix);

        % Method to return the majority center cone type for an RGC
        [theCenterConeTypeWeights, theCenterConeTypeNum, theMajorityConeType, ...
            theCenterConeTypes, theCenterConeIndices] = ...
            centerConeTypeWeights(obj, theRGCindex);

        % Method to return the L, M, and S net cone weights in the surround
        % If a non-empty theCenterConeIndices is passed the computed
        % surround cone type weights are computed after excluding the surround cones
        % that also feed into the RF center
        [theSurroundConeTypeWeights, theExclusiveSurroundConeTypeWeights, ...
         theSurroundConeTypes, theSurroundConeIndices] = ...
            surroundConeTypeWeights(obj, theRGCindex, theCenterConeIndices);

        % Method to return the index of the RGC best matching the target
        % criteria
        [targetCenterConesNumNotMatched, theCurrentRGCindex] = indexOfRGCNearPosition(obj, ...
            targetRGCposition, targetCenterConesNum, targetCenterConeMajorityType);

        % Method to return the indices of RGCs with a specific number of center cones
        indicesOfRGCsIdentified = indicesOfRGCsWithThisManyCenterCones(obj, targetCenterConesNum);

        % Method to eliminate RGCs with a specified number of center cones
        eliminateRGCsWithThisManyCenterConesNum(obj, targetCenterConesNum)

        % Stats on the center cones
        centerConesNumCases = centerConePoolingStats(obj);

        % Method to compute the spatial and chromatic variance cost for each RGC
        theCostComponentsMatrix = totalInputMaintenanceCost(obj, ...
            chromaticSpatialVarianceTradeoff, spatialVarianceMetric)

        % Method to overwrite the connectivity matrix and the RGCRF
        % positions (only used for debugging purposes)
        overwritePositionConnectivityDataForDebugPurposes(obj, theStageMetaDataStruct);

        % Setter for noiseFlag
        function set.noiseFlag(obj, val)
            assert((isempty(val))||ismember(val, mRGCMosaic.validNoiseFlags), ...
                sprintf('''%s'' is not a valid noise flag', val));
            obj.noiseFlag = val;
        end

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
        generateRFpositionsAndWireTheirCenters(obj,pResults);
        
        % Method to crop RGCs on the border
        cropRGCsOnTheBorder(obj);
    end % Private methods


    % Static methods
    methods (Static)
        % Method to render the cone pooling plot with a subregion based on
        % the cone indices and weights pooled by that subregion and return
        % the X,Y line weighting functions for that subregion
        [subregionLineWeightingFunctions, subregionContourData] = renderSubregionConePoolingPlot(ax, theConeMosaic, ...
            rgcRFposDegs, coneIndices, coneWeights, varargin);

        % Method to render the X,Y line weighting functions for a subregion
        renderSubregionConePoolingLineWeightingFunctions(ax, ...
            centerLineWeightingFunction, surroundLineWeightingFunction, ...
            sensitivityRange, horizontalAxisDirection, varargin);

        % Method to generate iso-sensitivity contour data for a 2D map at a specified set of levels
        cData = contourDataFromDensityMap(spatialSupportXY, zData, zLevels);

        % Method to list the available compute-ready mosaics
        availableComputeReadyMosaics(rgcMosaicType);

        % Method to load a compute-ready mosaic
        theComputeReadyMRGCmosaic = loadComputeReadyRGCMosaic(mosaicParams, opticsParams, retinalRFmodelParams);

        % Method to compute the extra degs to include the surrounds
        extraDegs = extraConeMosaicDegsForMidgetRGCSurrounds(eccentricityDegs, sizeDegs);

        % Method to generate  data for visualization of the PSF
        thePSFData = generateOpticsPSFdataForVisualization(theOI, visualizedWavelength, micronsPerDegree);

        % Method to generate RF center contours from the pooled cones when
        % the contourGenerationMethod is set to 'ellipseFitBasedOnLocalSpacing'
        contourData = subregionOutlineContourFromSpacing(...
            rgcPos, rgcRFradius,  ...
            xSupport, ySupport, spatialSupportSamples);

        % Method to generate RF center contours from the pooled cones when
        % the contourGenerationMethod is set to 'contourOfPooledConeApertureImage'
        [contourData, RFmap2D] = subregionOutlineContourFromPooledCones(...
            conePos, coneRc, poolingWeights, ...
            xSupport, ySupport, spatialSupportSamples)

        % Method to generate RF center contours from the pooled cones when
        % the contourGenerationMethod is set to 'ellipseFitToPooledConeApertureImage'
        contourData = subregionEllipseFromPooledCones(...
            conePos, coneRc, poolingWeights, ...
            xSupport, ySupport, spatialSupportSamples, centerSubregionContourSamples);

        % Method to fit an ellipse at some RFmap at the normalized level zLevel
        contourData = ellipseContourFromSubregionRFmap(xSupport, ySupport, RFmap, ...
            zLevel, centerSubregionContourSamples);
    end % Static methods

end



 