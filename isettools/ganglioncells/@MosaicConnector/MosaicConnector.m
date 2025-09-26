classdef MosaicConnector < handle
% Abstract class for all mosaic-connector-type objects:
% i.e., @mRGCConnector. All subclasses of
% MosaicConnector must implement the abstract methods defined here, 
% but the implementation may be different for different connectors.
% Also all subclasses  must specify the properties defined here.
% The data contained in these properties may differ from one
% class to another but their names must be the same.
% By using an abstract superclass we  enfore consistency amongst all 
% subclasses and also collect all common methods here. 
% 
%
% 8/07/2022  npc   Wrote it.
%
    % Public access.
    properties
        % Verbosity level (1 = minimal, 10 = max)
        verbosity

    end % Public properties

    % Constant properties
    properties (Constant)
        maxNeighborsNum = 6;                                    % max number of neighboring destination RF
        maxFractionOfRFspacingDuringInitialNeighborSearch = 0.6;
        maxNeighborNormDistance = 1.35;                         % max distance to search for neighboring destination RFs
        localSpacingFromCurrentCentroids = false;
    end

    % Protected properties. All @MosaicConnector subclasses can read these, 
    % but they cannot set them. 
    properties (SetAccess = protected)

        % The input source RF lattice
        sourceLattice;
        
        % The input destination RF lattice
        destinationLattice;
        
        % Indices of source RFs that are allowed to connect to a destination RF
        connectableSourceRFindices;

        % Wiring params struct
        customWiringParams;

        % Compute struct for computing local source-to-destination density ratios
        sourceToDestinationDensityRatioComputeStruct;

        % Centroids of destination RFs based on their current source RF inputs
        destinationRFcentroidsFromInputs;

        % Spacings of destination RFs based on their current centroids
        % (i.e., destinationRFcentroidsFromInputs), which are in turn based on their source RF inputs
        destinationRFspacingsFromCentroids = [];

        % Center of the source lattice
        sourceLatticeCenter = [];

        % Sparse [sourceRFsNum x destinationRFsNum] sparse  connectivity matrix 
        % To find which source RFs are connected to a targetDestinationRF:
        %  connectivityVector = full(squeeze(obj.connectivityMatrix(:, targetDestinationRF)));
        %  inputSourceRFIDs = find(connectivityVector > 0.01);
        connectivityMatrix = [];

        % Flag indicating whether to save meta data for each connectivity stage
        saveIntermediateConnectivityStagesMetaData;

        % Flag indicating whether to visualize connectivity as each connectivity stage
        visualizeConnectivityAtIntermediateStages;

        % Flags indicating whether to smooth the spacing noise( due to local positional jitter)
        % of the source and destination lattices
        smoothSourceLatticeSpacings;
        smoothDestinationLatticeSpacings;

        % Struct with wiring params
        wiringParams;

        % Cell array of meta data structs at each intermediate connectivity stage
        intermediateMetaDataStructs;

        % Cell array of figure handles for each intermediate connectivity stage
        intermediateFigureHandles = {};
    end % Write-protected 

    % The MosaicConnector subclass has no need to access these properties 
    % directly so they are protected
    properties (SetAccess = protected, GetAccess = protected)
    end % Fully protected


    % Abstract, public methods. Each subclass *must* implenent its own
    % version of all functions listed as abstract. If it does not, 
    % it cannot instantiate any objects.
    methods(Abstract)
        % Subclass-specific method to crop the source lattice
        % (depending on the destination lattice)
        cropSourceLattice(obj);

        % Subclass-specific method to crop the destination lattice
        % (depending on the source lattice). For example, when connecting
        % a cone mosaic (source) to an RGC mosaic (destination), we need
        % to crop the destination mosaic (RGC) so that cones extend
        % beyoind the edges of the RGC mosaic, to allow for surround inputs
        % to those RGCs at the edges of the RGC mosaic
        cropDestinationLattice(obj);

        % Subclass-secific method for returning the names of the different cost components
        costComponentNames = costComponentNames(obj);

        % Subclass-specific method to visualize the source lattice RFs
        visualizeSourceLatticeRFs(obj);

        % Subclass-specific method to visualize the statistic of the different cost components
        visualizeCostComponentStatistics(obj, ax1, ax2, theCostComponentsMatrix);

    end % Abstract methods

    % Public methods
    methods
        % Constructor
        function obj = MosaicConnector(sourceLattice, destinationLattice, varargin)
            % Parse input
            p = inputParser;
            p.addParameter('verbosity', 1);
            p.addParameter('densityRatioMapSamplingIntervalMicrons', 3, @isscalar);
            p.addParameter('connectableSourceRFindices', []);
            p.addParameter('visualizeConnectivityAtIntermediateStages', false, @islogical);
            p.addParameter('saveIntermediateConnectivityStagesMetaData', false, @islogical);
            p.addParameter('smoothSourceLatticeSpacings', true, @islogical);
            p.addParameter('smoothDestinationLatticeSpacings', true, @islogical);
            p.addParameter('wiringParams', [], @(x)(isempty(x)||isstruct(x)));
            
            % Execute the parser
            p.parse(varargin{:});
            obj.connectableSourceRFindices = p.Results.connectableSourceRFindices;
            obj.verbosity = p.Results.verbosity;
            obj.visualizeConnectivityAtIntermediateStages = p.Results.visualizeConnectivityAtIntermediateStages;
            obj.saveIntermediateConnectivityStagesMetaData = p.Results.saveIntermediateConnectivityStagesMetaData;
            obj.smoothSourceLatticeSpacings = p.Results.smoothSourceLatticeSpacings;
            obj.smoothDestinationLatticeSpacings = p.Results.smoothDestinationLatticeSpacings;
            obj.wiringParams = p.Results.wiringParams;

            % Validate source and destination lattices
            obj.validateInputLattice(sourceLattice, 'source');
            obj.validateInputLattice(destinationLattice, 'destination');

            % Generate the sourceToDestinationDensityRatioComputeStruct
            obj.generateSourceToDestinationDensityRatioComputeStruct(p.Results.densityRatioMapSamplingIntervalMicrons);

            % Crop the destination lattice - subclass specific
            obj.cropDestinationLattice();

            % Crop the source lattice - subclass specific
            obj.cropSourceLattice();
        end % Constructor

        % Method to update the intermediate meta data structs;
        updateIntermediateMetaDataStructs(obj, phaseDescriptor, netCostSequences, netTransfers);

        % Visualization methods
        hFig = visualizeCurrentConnectivity(obj, figNo, varargin);
        [hFig, ax, XLims, YLims] = visualizeInputLattices(obj, varargin);
        
        function set.sourceLattice(obj, s)
            obj.sourceLattice = s;
        end

        function set.destinationLattice(obj, s)
            obj.destinationLattice = s;
        end

        % Setter method for property verbosity
        function set.verbosity(obj, new_verbosity)
            obj.verbosity = new_verbosity;
            %fprintf('\nNew verbosity level: %d\n', obj.verbosity);
        end

        % Getter method for property verbosity
        function val = get.verbosity(obj)
            val = obj.verbosity;
        end

    end % public methods

    % Methods implemented in @MosaicConnector, which may be called its
    % subclasses, but are otherwise private
    methods (Access = protected)

        % Multi-step connection method
        connect(obj, varargin);

        % Compute the cost to maintain the current inputs for all destination RFs
        theCostComponentsMatrix = totalPoolingCosts(obj);

        % Update the destination RF centroids based on their inputs
        updateDestinationCentroidsFromInputs(obj, destinationRFList);

        % Update the destination RF spacings based on their centroids
        updateDestinationRFspacingsBasedOnCentroids(obj);
    end
    
    methods (Access = private)
        % Method to generate the sourceToDestinationDensityRatioComputeStruct
        generateSourceToDestinationDensityRatioComputeStruct(obj, samplingIntervalMicrons);

        % Method to compute the source:destination density ratio map
        densityRatioMap = sourceToDestinationDensityRatioMap(obj);

        % Stage 1 connection method
        connectSourceRFsToDestinationRFsBasedOnLocalDensities(obj);

        % Stage 2 connection method
        connectUnconnectedSourceRFsToClosestDestinationRF(obj);

        % Stage 3 connection method
        transferSourceRFsToZeroInputDestinationRFs(obj, varargin);
        removeZeroInputDestinationRFs(obj, indicesOfZeroInputDestinationRFs);

        % Stage 4/5 connection method
        theOptimizedDestinationRFindices = transferSourceRFsBetweenNearbyDestinationRFs(obj, varargin);

        % Method to compute sorted indices of destiation RFs based on their eccentricity &
        % optimization center
        [sortedIndices, sortedEccs] = sortDestinationRFsBasedOnOptimizationCenter(obj,unsortedIndices, varargin);

        % Input lattice validation method
        validateInputLattice(obj, theLattice, latticeName);
    end

    % Static methods
    methods (Static)
        [f,v] = facesAndVertices(positions, spacings, shapeOutline);

        [D,idx] = pdist2(A, B, varargin);
        convergenceIsAchieved = convergenceAchieved(netTotalCostSequence);
 
        visualizeConvergenceSequence(currentPass, ...
            costsMatrix, costsNames, ...
            netTransfers, maxPassesNum, plotTitle, figNo);

        % Method to determine the max number of cones that can be
        % transferred or swapped between nearby RF centers, depending on
        % the mosaic eccentricity
        [maxConesNumTransferred, maxConesNumSwapped] = coneTransferAndSwapLimits(mosaicEccDegs);

    end

end