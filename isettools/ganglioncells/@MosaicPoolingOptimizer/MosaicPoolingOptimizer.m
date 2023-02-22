classdef MosaicPoolingOptimizer < handle

    % Read-only properties 
    properties (SetAccess = private)
        % The RGC mosaic to be optimized
        theRGCMosaic;

        % The retinal RF model pooling params
        retinalRFmodelParams;

        % The nominal spatial sampling grid
        nominalSpatialSamplingGrid;

        % The sampling grids
        conesNumPooledByTheRFcenterGrid;
        targetRGCindicesWithLconeMajorityCenter;
        targetRGCindicesWithMconeMajorityCenter;
        visualSTFSurroundToCenterRcRatioGrid;
        visualSTFSurroundToCenterIntegratedSensitivityRatioGrid;
        
    end % Read only properties

    % Dependent properties
    properties (Dependent)
        % Number of grid nodes
        gridNodesNum;
    end


    properties (Constant)
        defaultRetinalRFmodelParams = struct(...
            'retinalConePoolingModel', 'arbitraryCenterConeWeights_doubleExpH1cellIndex1SurroundWeights', ...
            'H1cellIndex', 1, ...                                          
            'centerConnectableConeTypes', [cMosaic.LCONE_ID cMosaic.MCONE_ID], ...   % cone types that can connect to the RF center
            'surroundConnectableConeTypes', [cMosaic.LCONE_ID cMosaic.MCONE_ID] ... % cone types that can connect to the RF surround
        );

    end % Constants

    % Public methods (class interface)
    methods
        % Constructor
        function obj = MosaicPoolingOptimizer(theRGCmosaic, ...
            varargin)

            % Parse input
            p = inputParser;
            p.addRequired('theRGCmosaic', @(x)(isa(x, 'mRGCMosaic')));
            p.addParameter('minSpatialSamplingDegs', 0.25, @isnumeric);
            p.addParameter('samplingScheme', 'hexagonal', @(x)(ismember(x, {'hexagonal', 'rectangular'})));
            p.parse(theRGCmosaic, varargin{:});

            obj.theRGCMosaic = p.Results.theRGCmosaic;
            obj.retinalRFmodelParams = obj.defaultRetinalRFmodelParams;

            % Generate the nominal multifocal spatial sampling grid
            obj.generateNominalSpatialSamplingGrid(p.Results.samplingScheme);

            % Generate the full multifocal sampling grids
            obj.generateSamplingGrids(p.Results.minSpatialSamplingDegs);

            obj.visualizeSamplingGrids();

        end % Constructor

        % Visualize Sampling grids
        visualizeSamplingGrids(obj);

        % Fit the optimizer
        fit(obj, varargin);

        % Compute cone mosaicSTFs resources for the fitter at a single grid node
        computeConeMosaicSTFresponses(obj, gridNode, stimSizeDegs, responsesFileName, varargin);

        % Getter for dependent property gridNodesNum
        function val = get.gridNodesNum(obj)
            val = numel(obj.targetRGCindicesWithLconeMajorityCenter);
        end


    end % Public methods

    % Private methods
    methods (Access=private)
        % Method to generate the nominal spatial sampling grid
        generateNominalSpatialSamplingGrid(obj, samplingScheme);

        % Method to generate the full sampling grids
        generateSamplingGrids(obj, minSpatialSamplingDegs);

        % Method to setup the parameters and the display for conducting an
        % STF mapping experiment.
        [stimParams, thePresentationDisplay] = setupSTFmappingExperiment(obj, ...
            sceneFOVdegs, retinalImageResolutionDegs)

        % Method to return the center majority cone types
        [theCenterConeTypeWeights, theCenterConeTypeNum, theMajorityConeType] = centerConeTypeWeights(obj, theRGCindex);
    end

    % Static methods
    methods (Static)
        dropboxDir = localDropboxPath();
    end % Static methods

end
