% Subclass of @MosaicConnector, for connecting mRGClattice to a cone lattice
%

classdef coneToMidgetRGCConnector < MosaicConnector

    % Public properties (specific to the coneToMidgetRGCConnector) 
    properties
        
    end
    
    % Constant properties
    properties (Constant)

        defaultWiringParams = struct(...
            'optimizationCenter', 'latticeCenter', ...            % {'latticeCenter', 'origin'}
            'chromaticSpatialVarianceTradeoff', 1.0, ...          % [0: minimize chromatic variance 1: minimize spatial variance]
            'rfCentroidOverlapPenaltyFactor', 1, ...              % Penalty for overlapping centroids
            'destinationRFoverlapRatio', 0.0, ...                 % overlap of midgetRGCRFs (0 = no overlap)
            'spatialVarianceMetric', 'spatial variance', ...      % choose between {'maximal interinput distance', 'spatial variance'}
            'maxMeanConeInputsPerRGCToConsiderSwapping', 7, ...   % Do cone swapping only if the mean cones/RGC less than or equal to this number
            'maxNumberOfConesToSwap', 4, ...                      % Only swap up to this many cones
            'maxPassesNum', 15 ...     
        );

    end

    % --- PRIVATE PROPERTIES ----------------------------------------------
    properties (Access = private)            
        coneTypeInfoIsAvailable = false;
    end
    % --- END OF PRIVATE PROPERTIES ---------------------------------------
    
    
    % Public methods
    methods
        % Constructor
        function obj = coneToMidgetRGCConnector(...
                sourceLattice, destinationLattice, varargin) 

            p = inputParser;
            p.addParameter('verbosity', 1);
            p.addParameter('generateProgressVideo', false, @islogical);
            p.addParameter('coneIndicesToBeConnected', []);
            p.addParameter('visualizeConnectivityAtIntermediateStages', false, @islogical);
            p.addParameter('saveIntermediateConnectivityStagesMetaData', false, @islogical);
            p.addParameter('smoothSourceLatticeSpacings', true, @islogical);
            p.addParameter('smoothDestinationLatticeSpacings', true, @islogical);

            p.addParameter('maxNeighborNormDistance', MosaicConnector.maxNeighborNormDistance, @isscalar);
            p.addParameter('maxNeighborsNum', MosaicConnector.maxNeighborsNum, @isscalar);
            p.addParameter('maxConeInputsPerRGCToConsiderTransferToNearbyRGCs', MosaicConnector.maxSourceInputsToConsiderTransferToNearbyDestinationRF, @(x)(isscalar(x)&&(x>=MosaicConnector.minSourceInputsToConsiderTransferToNearbyDestinationRF)));
            p.addParameter('maxMeanConeInputsPerRGCToConsiderSwappingWithNearbyRGCs',MosaicConnector.maxMeanSourceInputsToConsiderSwappingWithNearbyDestinationRF, @(x)(isscalar(x)&&(x>1)));
            p.addParameter('chromaticSpatialVarianceTradeoff', coneToMidgetRGCConnector.defaultWiringParams.chromaticSpatialVarianceTradeoff, @(x)(isscalar(x)&&(x>=0)&&(x<=1)));  % [0: minimize chromatic variance 1: minimize spatial variance]
            p.addParameter('optimizationCenter', coneToMidgetRGCConnector.defaultWiringParams.optimizationCenter, @(x)(ismember(x, {'patchCenter', 'origin'})));

            % Execute the parser
            p.parse(varargin{:});
            
            customWiringParams = coneToMidgetRGCConnector.defaultWiringParams;
            customWiringParams.maxNeighborNormDistance = p.Results.maxNeighborNormDistance;
            customWiringParams.maxNeighborsNum = p.Results.maxNeighborsNum;

            customWiringParams.maxConeInputsPerRGCToConsiderTransferToNearbyRGCs = p.Results.maxConeInputsPerRGCToConsiderTransferToNearbyRGCs;
            customWiringParams.maxMeanConeInputsPerRGCToConsiderSwapping = p.Results.maxMeanConeInputsPerRGCToConsiderSwappingWithNearbyRGCs;

            customWiringParams.chromaticSpatialVarianceTradeoff = p.Results.chromaticSpatialVarianceTradeoff;
            customWiringParams.optimizationCenter = p.Results.optimizationCenter;
           
            generateProgressVideo = p.Results.generateProgressVideo;

            % Initialize the super-class constructor.
            obj = obj@MosaicConnector(...
                sourceLattice, ...
                destinationLattice, ...
                'verbosity', p.Results.verbosity, ...
                'connectableSourceRFindices', p.Results.coneIndicesToBeConnected, ...
                'saveIntermediateConnectivityStagesMetaData', p.Results.saveIntermediateConnectivityStagesMetaData, ...
                'visualizeConnectivityAtIntermediateStages', p.Results.visualizeConnectivityAtIntermediateStages, ...
                'smoothSourceLatticeSpacings', p.Results.smoothSourceLatticeSpacings, ...
                'smoothDestinationLatticeSpacings', p.Results.smoothDestinationLatticeSpacings, ...
                'wiringParams', customWiringParams);

            % Update our own variables
            if (isfield(sourceLattice, 'metaData')) && ...
                (isfield(sourceLattice.metaData, 'coneTypes')) && ...
                (isfield(sourceLattice.metaData, 'coneTypeIDs'))
                    obj.coneTypeInfoIsAvailable = true;
            end

            % Connect the 2 mosaics
            obj.connect('generateProgressVideo', generateProgressVideo);
        end % Constructor


    end % Public methods

    % Implementations of required -- Public -- Abstract methods defined in the MosaicConnector interface   
    methods
        % coneToMidgetRGCconnector - specific method to crop the source lattice
        % (depending on the destination lattice)
        cropSourceLattice(obj);

        % coneToMidgetRGCconnector - specific-specific method to crop the destination lattice
        % (depending on the source lattice)
        cropDestinationLattice(obj);

        % coneToMidgetRGCconnector -specific method to compute the cost
        % components to maintain a set of inputs (cones)
        theCostComponents = inputMaintenanceCost(obj, inputIndices, inputWeights, destinationRFspacing);

        % Return names of the different cost components
        costComponentNames = costComponentNames(obj);

        

        % Subclass-secific method for computing the various cost components
        % to maintain the overlap between two RGCs
        theCostComponents = overlappingDestinationRFCost(obj, ...
            theRGCindex, ...
            theRGCinputIndices, theRGCinputWeights, ...
            theNearbyRGCindex, ...
            theNearbyRGCinputIndices, theNearbyRGCinputWeights ...
            );

        % coneToMidgetRGCconnector - specific method to visualize the
        % source lattice RFs (i.e., the cone RFs)
        visualizeSourceLatticeRFs(obj, ax, coneOutline, varargin);

        % coneToMidgetRGCconnector - specific method to visualize the
        % statistic of the different cost components
        visualizeCostComponentStatistics(obj, ax1, ax2, theCostComponentsMatrix);

    end % Implementations of required -- Public -- Abstract methods defined in the MosaicConnector interface
    

    methods (Access = private) 

    end % Private methods 

    methods (Static)
        cost = chromaticVarianceCost(inputWeights, inputConeTypes);
        cost = spatialVarianceCost(spatialVarianceMetric, inputWeights, inputPositions, destinationRFspacing);
    end % Static methods

end


