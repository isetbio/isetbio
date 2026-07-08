% Subclass of @MosaicConnector, for connecting mRGClattice to a cone lattice
%

classdef coneToMidgetRGCConnector < MosaicConnector

    % Public properties (specific to the coneToMidgetRGCConnector) 
    properties
        
    end
    
    % Constant properties
    properties (Constant)

        defaultWiringParams = struct(...
            'optimizationCenter', 'latticeCenter', ...            % {'latticeCenter', 'origin', or, 'localSpacingBased'}
            'spatialChromaticUniformityTradeoff', 1.0, ...        % [0: maximize RF center chromatic uniformity, 1: maximize RF center spatial uniformity]
            'destinationRFoverlapRatio', 0.0, ...                 % overlap of midgetRGCRFs (0 = no overlap)
            'maxSourceInputsToConsiderTransferToNearbyDestinationRF', 30, ...
            'maxConeInputsPerRGCToConsiderSwappingWithNearbyRGCs', 4, ...                      % Only swap cones for  RGCs with up to this many input cones
            'maxPassesNum', 10 ...     
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
            p.addParameter('localSpacingFromCurrentCentroids', MosaicConnector.localSpacingFromCurrentCentroids);
            p.addParameter('spatialChromaticUniformityTradeoff', coneToMidgetRGCConnector.defaultWiringParams.spatialChromaticUniformityTradeoff, @(x)(isscalar(x)&&(x>=0)&&(x<=1)));  % [0: maximize RF center chromatic uniformity, 1: maximize RF center spatial uniformity]
            p.addParameter('optimizationCenter', coneToMidgetRGCConnector.defaultWiringParams.optimizationCenter, @(x)(ismember(x, {'patchCenter', 'origin', 'localSpacing', 'localConeToRGCdensityRatio'})));
            
            p.addParameter('maxConeInputsPerRGCToConsiderTransferToNearbyRGCs', coneToMidgetRGCConnector.defaultWiringParams.maxSourceInputsToConsiderTransferToNearbyDestinationRF, @isscalar);
            p.addParameter('maxConeInputsPerRGCToConsiderSwappingWithNearbyRGCs', coneToMidgetRGCConnector.defaultWiringParams.maxConeInputsPerRGCToConsiderSwappingWithNearbyRGCs, @isscalar);
            p.addParameter('maxPassesNum', coneToMidgetRGCConnector.defaultWiringParams.maxPassesNum, @(x)(isscalar(x)&&(x>=1)));

            % Execute the parser
            p.parse(varargin{:});
            
            customWiringParams = coneToMidgetRGCConnector.defaultWiringParams;
            customWiringParams.maxNeighborNormDistance = p.Results.maxNeighborNormDistance;
            customWiringParams.maxNeighborsNum = p.Results.maxNeighborsNum;
            customWiringParams.maxConeInputsPerRGCToConsiderTransferToNearbyRGCs = p.Results.maxConeInputsPerRGCToConsiderTransferToNearbyRGCs;
            customWiringParams.maxConeInputsPerRGCToConsiderSwappingWithNearbyRGCs = p.Results.maxConeInputsPerRGCToConsiderSwappingWithNearbyRGCs;
            customWiringParams.localSpacingFromCurrentCentroids = p.Results.localSpacingFromCurrentCentroids;
            customWiringParams.spatialChromaticUniformityTradeoff = p.Results.spatialChromaticUniformityTradeoff;
            customWiringParams.optimizationCenter = p.Results.optimizationCenter;
            customWiringParams.maxPassesNum = p.Results.maxPassesNum;

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

        % Return names of the different cost components
        costComponentNames = costComponentNames(obj);

        % coneToMidgetRGCconnector - specific method to visualize the
        % source lattice RFs (i.e., the cone RFs)
        visualizeSourceLatticeRFs(obj, ax, coneOutline, varargin);

        % coneToMidgetRGCconnector - specific method to visualize the
        % statistic of the different cost components
        visualizeCostComponentStatistics(obj, ax1, ax2, theCostComponentsMatrix);

    end % Implementations of required -- Public -- Abstract methods defined in the MosaicConnector interface
    

    methods (Access = protected) 
        theTotalPoolingCosts = totalPoolingCosts(obj);
    end % Private methods 

    methods (Static)
        cost = spectralUniformityCost(inputConeTypes, inputConeWeights);
        [theCentroidSpacingCost, theCenterConeNumerosityDifferential, theCentroidOverlapCost, theVarianceCost] = spatialCompactnessCost(...
            theTargetRFSourceRFpositions, theDonorRFSourceRFpositions, ...
            theTargetRFSourceRFconnectionWeights, theDonorRFSourceRFconnectionWeights, ...
            theTargetDestinationRFspacing, theDonorDestinationRFspacing);
    end % Static methods

end


