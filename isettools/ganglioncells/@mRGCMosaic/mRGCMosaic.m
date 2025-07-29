classdef mRGCMosaic < handle
% Create an mRGCMosaic mosaic object

    % Public properties
    properties  (GetAccess=public, SetAccess=public)
        % Name of the mosaic
        name = 'unnamed MRGC mosaic';

        % noise flag
        noiseFlag = 'random';

        % random seed
        randomSeed = 1;

        % vMembrane noise sigma
        vMembraneGaussianNoiseSigma = 0.07;

        % Scaling factor for temporal equivalent eccentricity
        temporalEquivantEccentricityScalingMode = 'ISETBioMosaicsBased';
    end % Public properties


    % Constant properties
    properties (Constant)
        % Valid noise flags
        validNoiseFlags = {'none', 'frozen', 'random'};
        defaultRandomSeed = 123;
        
        validOpticsModifications = {'adaptiveOptics6MM', 'adaptiveOptics6MMwithLCA', 'diffractionLimited', 'centralRefractionCorrected', 'maxStrehlRatio', 'customRefraction', 'nativeOpticsNoStrehlRatioOptimization'};

        % max surround support factor (used during optimizing surround)
        maxSurroundSupportFactor = 7.0;   % max surround radius: this factor x C&K surroundCharacteristicRadius

        % Min surround weight for inclusion in response computation
        minSurroundWeightForInclusionInComputing = 1e-4;


        % default minConeWeight when none is specified.
        % This is the amplitude of a Gaussian at 1 sigma, the weight at which nearby
        % RGCs abut according to Gauthier et al (2008). 
        % "Uniform Signal Redundancy of Parasol and Midget Ganglion Cells in Primate Retina"
        % Since the abut at 1 sigma, the sensitivity should be exp(-1/2). 
        sensitivityAtPointOfOverlapFromGauthierRFmappingExperiment = exp(-0.5);

        % However, since we are doing the overlap with no surrounds, and the in-vitro measurements are
        % with the surrounds (obviously), this sensitivity should be amplified by some factor.
        % Different attempts:
        %    0.85 -> way to much overlap and RF centers quite large (first attempt)
        %    0.55 -> a little too little overlap; 
        %    0.70 -> still too much overlap (the .1 eccentricity mosaics were done with 0.7); 
        %    0.60 -> probably the best value, going with that
        amplificationInCenterOnlySensitivityCausedByInactiveSurrounds = 0.6; 

        % Sensitivity for RF center only operations
        sensitivityAtPointOfOverlap = mRGCMosaic.sensitivityAtPointOfOverlapFromGauthierRFmappingExperiment * ...
                                      mRGCMosaic.amplificationInCenterOnlySensitivityCausedByInactiveSurrounds;

        % Other RF center overlap parameters:
        % The max exponent, 10 (flat-top). The min is 2 (regular Gaussian)
        superGaussianMaxExponent = 2 + 8;
        % The ecc-dependent variation in superGaussian exponent: steepness
        superGaussianExponentSigmoidSteepnessExponent = 10; 
        % The ecc-dependent variation in superGaussian exponent: 'c50' point
        superGaussianExponentSigmoidEccDegs50 = 8;

        % 1% min sensitivity for inclusion of divergence cone connections
        minSensitivityForInclusionOfDivergentConeConnections = 1/100;

        % 10% which is the noise floor in the measurements of cone weights according to RF center overlap according to Greg Field
        % and which is what they used in Fig 4 of their 2010 paper.
        minRFcenterConeWeightIncludedToMatchFigure4OfFieldEtAl2010 = 10/100;

        % In the supplementary section of Field et al (2010) paper it is stated that all cells were recorded at
	    % an eccentricity of 6.75 mm along the horizontal meridian.
	    % 6.75 mm in macaque retina corresponds to 25 degs in human retina
        eccentricityMMsForFigure4OfFieldEtAl2010 = 6.75;

        % Convert to degs in macaque retina
	    eccentricityDegsMacaqueRetinaForFigure4OfFieldEtAl2010 = ...
            RGCMosaicConstructor.helper.convert.eccentricityInMacaqueRetina(...
               'MMsToDegs', mRGCMosaic.eccentricityMMsForFigure4OfFieldEtAl2010);

	    % Convert to degs in human retina
	    eccentricityDegsHumanRetinaForFigure4OfFieldEtAl2010 = ...
            RGCMosaicConstructor.helper.convert.eccentricityDegsBetweenSpecies(...
            'MacaqueRetinaToHumanRetina', RGCMosaicConstructor.helper.convert.eccentricityInMacaqueRetina(...
                                               'MMsToDegs', mRGCMosaic.eccentricityMMsForFigure4OfFieldEtAl2010));


        % Valid methors for RF subregion contour generation
        validRFsubregionContourGenerationMethods = {...
            'ellipseFitBasedOnLocalSpacing', ...
            'contourOfPooledConeApertureImage', ...
            'ellipseFitToPooledConeApertureImage', ...
            'ellipseFitToPooledConePositions'};

        validTemporalEquivalentScalingModes = {'ISETBioMosaicsBased', 'WatanabeRodieckBased'};
    end  % Constant properties


    % Read-only properties
    properties (GetAccess=public, SetAccess=private)

        % Eccentricity [x,y] of the center of the mRGC mosaic (degs)
        eccentricityDegs;

        % Eccentricity [x,y] of the center of the mRGC mosaic (microns)
        eccentricityMicrons;

        % Size of the mRGC mosaic [width, height]
        sizeDegs;

        % Either 'right eye' or 'left eye'
        whichEye;

        % The params that were used to determine the center-connectivity
        rfCenterConnectivityParams;

        % The input cone mosaic (@cMosaic)
        inputConeMosaic;

        % The number of RGCs
        rgcsNum;

        % [n x 2] matrix of rgc rf positions, in microns & degrees
        rgcRFpositionsMicrons;
        rgcRFpositionsDegs;

        % [n x 1] vector of rgc rf spacings, in microns & degrees
        rgcRFspacingsMicrons;
        rgcRFspacingsDegs;

        % Positions and spacings of the source (not connected to cones mRGC lattice)
        rgcRFpositionsDegsOfSourceLattice;
        rgcRFspacingsDegsOfSourceLattice;

        % Flag indicating whether RF centers have been wired to overlap
        rgcRFcentersOverlap = false;

        % Struct with rfCenter overlapping params
        rfCenterOverlapParams = struct();


        % Computed params
        % Sparse [conesNum x rgcRFsNum] connectivity matrix for the centers
        % To find which cones are connected to an rgcRF:
        %  connectivityVector = full(squeeze(rgcRFcenterConeConnectivityMatrix(:, rgcRF)));
        %  inputConeIndices = find(abs(connectivityVector) > 0.0001);
        rgcRFcenterConeConnectivityMatrix = [];
        rgcRFcenterConeConnectivityMatrixNoOverlap = [];
        exclusivelyConnectedInputConeIndicesNum = [];

        % Sparse [conesNum x rgcRFsNum] connectivity matrix for the surrounds
        rgcRFsurroundConeConnectivityMatrix = [];

        % The response gain for each RGC
        responseGains = [];

        % Struct with all the params that were used to generate the rgcRFsurroundConeConnectivityMatrix
        rfSurroundConnectivityParams = [];

        % Struct with variance params for the surround
        surroundVarianceInComputeReadyMosaic = [];

        visualizationCache = [];
    end %  Read-only properties

    % Dependent properties
    properties (Dependent)
    end % Dependent properties


    % Private properties
    properties (GetAccess=private, SetAccess=private)
    end % private properties

    % Public methods
    methods
        % Constructor
        function obj = mRGCMosaic(varargin)

            % If the input arguments are in the form of a struct, we
            % leave it alone.
            %
            % The key/val parameters are all defined as lower case
            % with no spaces. The ieParamFormat function allows the
            % user to specify parameters that include spaces and mixed
            % case.  It forces to lower case and removes space.  This
            % improves readability.
            if ~isempty(varargin) && ~isstruct(varargin{1})
                varargin = ieParamFormat(varargin);
            end

            % Parse input
            p = inputParser;
            p.addParameter('withcomponents', [], @(x)(isempty(x) || (isstruct(x))));
            p.addParameter('eccentricitydegs', [], @(x)(isnumeric(x) && (numel(x) == 2)));
            p.addParameter('positiondegs', [], @(x)(isnumeric(x) && (numel(x) == 2)));
            p.addParameter('sizedegs', [0.4 0.4], @(x)(isnumeric(x) && (numel(x) == 2)));
            p.addParameter('whicheye', 'right eye', @(x)(ischar(x) && (ismember(x, {'left eye', 'right eye'}))));
            p.parse(varargin{:});

            if (~isempty(p.Results.withcomponents))
                obj.initializeWithComponents(p.Results.withcomponents);
                fprintf('Initializing with passed components. All other key/value pair values are ignored.\n')
                return;
            end
        end % Constructor

        % Method to generate optics that are native to the MRGCmosaic
        [theOI, thePSF, theOptimalStrehlRatioDefocusDiopters, theOptimalStrehlRatio] = nativeOI(obj, varargin);

        % Method to compute  the spatiotemporal response of the mRGCMosaic given the response of its input cone  mosaic
        [noiseFreeMRGCresponses, noisyMRGCresponseInstances, responseTemporalSupportSeconds] = compute(obj, ...
            theInputConeMosaicResponse, theInputConeMosaicResponseTemporalSupportSeconds, varargin);

        % Method to generate noisy response instances
        noisyMRGCresponseInstances = noisyResponseInstances(obj, noiseFreeMRGCresponses, varargin);

        % Method to visualize the mRGCMosaic
        [hFig, ax] = visualize(obj, varargin);

        % Method to generate the visualization cache
        generateVisualizationCache(obj, xSupport, ySupport, ...
            centerSubregionContourSamples, contourGenerationMethod, ...
            visualizedRGCindices, minConedWeightIncluded, varargin);

        % Method to crop the mosaic
        theCroppedMosaic = cropToSizeAtEccentricity(obj, sizeDegs, eccentricityDegs, varargin);

        % Method to generate RF center overlap
        generateRFcenterOverlap(obj, rfCenterOverlapParamsStruct, varargin);

        % Method to set the  rgcRFsurroundConeConnectivityMatrix, the gains, the rfSurroundConnectivityParams
        bakeSurroundConeConnectivityMatrixAndFreeze(obj, rgcRFsurroundConeConnectivityMatrix, rfSurroundConnectivityParams, surroundVarianceInComputeReadyMosaic);

        % Method to compute the spatial and chromatic uniformity of cone inputs to the RF centers
        [theSpatialCompactnessCosts, theSpectralUniformityCosts, ...
        atLeastOneConeInputWithGreaterThanMinWeightIncluded, theInputConeNumerosityDifferentials, ...
        theCentroidSpacingCosts, theSpatialVarianceCosts] = ...
            rfCenterSpatioChromaticCosts(obj, varargin);

        % Method to compute the surround cone purities for all a set or all RGCs in the mosaic
        [surroundConePurities, centerConeDominances, centerConeNumerosities, centerConePurities] = surroundConePurities(obj, theRGCindices, surroundConeSelection);

        % Method to return cone connectivity stats for theSubregion ('center' or 'surround') oftheRGCindex
        s = singleCellConnectivityStats(obj, theRGCindex, theSubregion, varargin);

        % Method to return connectivity stats for the RF centers of all mRGCs within the mosaic
        [allRGCCenterConesNum, allRGCCenterDominantConeTypes, allRGCCenterRelativeConeWeights] = ...
            allRFcenterConnectivityStats(obj, varargin);

        % Method to return the indices of all RGCs with a specified range in the following properties:
        % - center cone numerocity 
        % - surround purity (for L-center RFs, surround purity is: M/(L+M), and correspondingly for M-center RFs)
        % - radial eccentricity
        % - center purity (for L-center RFs, center purity is: M/(L+M), and correspondingly for M-center RFs)

        [targetRGCindices, surroundConePurities, centerConeDominances, ...
          centerConeNumerosities, centerConePurities] = indicesOfRGCsWithinTargetedPropertyRanges(obj, ...
            targetedCenterConeNumerosityRange, ...
            targetedSurroundPurityRange, ...
            targetedRadialEccentricityRange, ...
            targetedCenterPurityRange)

        % Method to return the indices of all RGCs with a specificed number of center input cones
        theRGCIndices = indicesOfRGCsWithTargetCenterConesNumInRange(obj,targetCenterConesNumRange, varargin);

        % Method to return the indices of all RGCs within a rectangular region
        theRGCIndices = indicesOfRGCsWithinROI(obj, roiCenterDegs, roiSizeDegs);

        % Method to return the indices of all RGCs that are near targetRGC
        [nearbyRGCindices, distancesToNearbyRGCs] = indicesOfNearbyRGCs(obj, targetRGC, varargin);

        % Method to return the indices and distances of the maxNeighborsNum nearest RGCs to theTargetRGCindices
        % nearbyRGCindices: 
        %   [maxNeighborsNum x numel(theTargetRGCindices)] matrix of indices of nearby RGCs to theTargetRGCindices
        % distanceMatrix: 
        %   [maxNeighborsNum x numel(theTargetRGCindices)] matrix of distances of nearby RGCs to theTargetRGCindices
        [nearbyRGCindices, distancesToNearbyRGCs] = neighboringRGCsToTargetRGCs(obj, theTargetRGCindices, maxNeighborsNum)

        % Methods to go back and forth between angular and linear eccentricities
        angularEccDegs = linearEccMMsToAngularEccDegs(obj, linearEccMM);
        linearEccMMs = angularEccDegsToLinearEccMMs(obj, angularEccDegs);

        % Return temporal equivalent eccentricities (xy) for an array of xy positions
        temporalEquivalentEccDegs = temporalEquivalentEccentricityForEccXYDegs(obj, xyEccDegs);

        % Setters
        function set.temporalEquivantEccentricityScalingMode(obj, val)
            assert((ischar(val))&&ismember(val, mRGCMosaic.validTemporalEquivalentScalingModes), ...
                sprintf('''%s'' is not a valid temporalEquivantEccentricityScalingMode', val));
            obj.temporalEquivantEccentricityScalingMode = val;
        end

        function set.noiseFlag(obj, val)
            assert((ischar(val))&&ismember(val, mRGCMosaic.validNoiseFlags), ...
                sprintf('''%s'' is not a valid noise flag', val));
            obj.noiseFlag = val;
        end

        function set.randomSeed(obj, val)
            assert(isnumeric(val), sprintf('randomSeed must be a numeral'));
            obj.randomSeed = val;
        end

        function set.vMembraneGaussianNoiseSigma(obj, val)
            assert(isnumeric(val)&&(val>=0), sprintf('rvMembraneGaussianNoiseSigma must be a positive numeral'));
            obj.vMembraneGaussianNoiseSigma = val;
        end

    end % Public methods


    % Private methods
    methods (Access = private)

        % Method to initialize an mRGCMosaic with previously-computed components
        % These components are generated by the RGCMosaicConstructor
        function initializeWithComponents(obj,pStruct)
            propertyNames = fieldnames(pStruct);
            for iProperty = 1:numel(propertyNames)
                thePropertyName = propertyNames{iProperty};
                if (isprop(obj,thePropertyName))
                    fprintf('Generating mRGCMosaic (property: ''%s'').\n', thePropertyName);
                    obj.(thePropertyName) = pStruct.(thePropertyName);
                else
                    error('mRGCMosaic does not include a ''%s'' property.\n', thePropertyName);
                end
            end
        end
 
        % Method to render the cone pooling RF map of a single RGC. Called
        % by visualize(), using the following key-value pairs
        % 'visualizedRGCindices', [some list of RGCs], ...
        % 'singleRGCconePoolingRFmaps', true, 
        hFig = renderConePoolingRFmap(obj, theRGCindex, varargin);

    end % private methods


    % Static methods
    methods (Static)
        % Method to list available prebaked mRGCMosaics
        [themRGCmosaicFileNames, prebakedMRGCMosaicDir] = listPrebakedMosaics();

        % Method to load a prebaked mRGCMosaic and associated optics
        [theMRGCmosaic, theOI] = loadPrebakedMosaic(mosaicParams, opticsParams);

        % Method to compute extra support needed to account for surrounds at a given eccentricity
        extraDegs = extraSupportDegsForMidgetRGCSurrounds(eccentricityDegs, sizeDegs);

        % Method to render the input cone mosaic pooling map within a subregion (center or surround)
        subregionLineWeightingFunctions = renderInputConeMosaicSubregionPoolingMap(theAxes, ...
            theInputConeMosaic, theContourData, ...
            theSubregionConeIndices, theSubregionConeWeights, flatTopSaturationLevel, ...
            spatialSupportCenterDegs, spatialSupportTickSeparationArcMin, ...
            scaleBarDegs, gridless, noXLabel, noYLabel, varargin);

        % Method to compute the cone mosaic pooling map 
        [retinalSubregionConeMap, retinalSubregionConeMapFlatTop] = retinalSubregionConeMapFromPooledConeInputs(...
            theConeMosaic, theConeIndices, theConeWeights, spatialSupportDegs, flatTopSaturationLevel)

        % Method to fit an ellipse to a subregion defined by pooled cones
        contourData = subregionEllipseFromPooledCones(conePos, coneSpacings, poolingWeights, ...
                xSupport, ySupport, spatialSupportSamples, centerSubregionContourSamples);

        % Method to fit an ellipse to the pooled cones
        [contourData, theCenter, xAlpha, yAlpha, theRotationRadians] = subregionEllipseFromPooledConePositions(...
                         theConePositions, ellipseContourPoints, ellipseContourAngles, pLevel, maxNumberOfConesOutsideContour);
    end % Static methods
    
end
