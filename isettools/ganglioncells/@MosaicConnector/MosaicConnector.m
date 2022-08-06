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
% 3/27/2014  npc   Wrote it.
%
    % Public access.
    properties
        % Verbosity level (1 = minimal, 10 = max)
        verbosity

    end % Public properties

    % Protected properties. All @MosaicConnector subclasses can read these, 
    % but they cannot set them. 
    properties (SetAccess = protected)

        % The input source RF lattice
        sourceLattice;
        
        % The input destination RF lattice
        destinationLattice;
        
        % Compute struct for computing local source-to-destination density ratios
        sourceToDestinationDensityRatioComputeStruct;

        % Indices of source RFs that are allowed to connect to a destination RF
        connectableSourceRFindices;

        % Centroids of destination RFs based on the current source RF inputs
        destinationRFcentroidsFromInputs;

        % Sparse [sourceRFsNum x destinationRFsNum] sparse  connectivity matrix 
        % To find which source RFs are connected to a targetDestinationRF:
        %  connectivityVector = full(squeeze(obj.connectivityMatrix(:, targetDestinationRF)));
        %  inputSourceRFIDs = find(connectivityVector > 0.01);
        connectivityMatrix = [];


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
            % Execute the parser
            p.parse(varargin{:});

            obj.connectableSourceRFindices = p.Results.connectableSourceRFindices;
            obj.verbosity = p.Results.verbosity;

            % Assign and validate source and destination lattices
            obj.sourceLattice = sourceLattice;
            [obj.sourceLattice.RFspacingsMicrons, ...
             obj.sourceLattice.nearbyRFindices] = RGCmodels.Watson.convert.positionsToSpacings(obj.sourceLattice.RFpositionsMicrons);
            % Smooth source lattice spacings
            obj.sourceLattice.RFspacingsMicrons = MosaicConnector.smoothSpacings(...
                obj.sourceLattice.RFspacingsMicrons, ...
                obj.sourceLattice.nearbyRFindices);

            obj.destinationLattice = destinationLattice;
            [obj.destinationLattice.RFspacingsMicrons, ...
             obj.destinationLattice.nearbyRFindices] = RGCmodels.Watson.convert.positionsToSpacings(obj.destinationLattice.RFpositionsMicrons);
            % Smooth destination lattice spacings
            obj.destinationLattice.RFspacingsMicrons = MosaicConnector.smoothSpacings(...
                obj.destinationLattice.RFspacingsMicrons, ...
                obj.destinationLattice.nearbyRFindices);

            obj.validateInputLattice(obj.sourceLattice, 'source');
            obj.validateInputLattice(obj.destinationLattice, 'destination');

            % Generate the sourceToDestinationDensityRatioComputeStruct
            obj.generateSourceToDestinationDensityRatioComputeStruct(p.Results.densityRatioMapSamplingIntervalMicrons);

            % Crop the destination lattice - subclass specific
            obj.cropDestinationLattice();

            % Crop the source lattice - subclass specific
            obj.cropSourceLattice();

            % Step1. Connect mosaics based on the sourceToDestinationDensityRatio
            
            % Initialize centroids. No inputs so set them all to inf
            destinationRFsNum = size(obj.destinationLattice.RFpositionsMicrons,1);
            obj.destinationRFcentroidsFromInputs = inf(destinationRFsNum,2);

            obj.connectSourceRFsToDestinationRFsBasedOnLocalDensities();
            
        end % Constructor


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

    % The class-user has no need to call these methods, but our subclasses may.
    methods (Access = protected)
        
    end
    
    methods (Access = private)

        % Method to generate the sourceToDestinationDensityRatioComputeStruct
        generateSourceToDestinationDensityRatioComputeStruct(obj, samplingIntervalMicrons);

        % Method to compute the source:destination density ratio map
        densityRatioMap = sourceToDestinationDensityRatioMap(obj);

        % Stage1 connection methods
        connectSourceRFsToDestinationRFsBasedOnLocalDensities(obj);
        updateDestinationCentroidsFromInputs(obj, nearestDestinationRFIndices);

        % Input lattice validation method
        function validateInputLattice(~, theLattice, latticeName)
            % Must be a struct
            assert(isstruct(theLattice), 'The %s lattice must be a struct.', latticeName);
            
            % Must have fields
            assert(isfield(theLattice, 'DegsToMMsConversionFunction')&&(isa(theLattice.DegsToMMsConversionFunction,'function_handle')), 'The %s lattice must have a ''DegsToMMsConversionFunction'' function handle field.', latticeName);
            assert(isfield(theLattice, 'MMsToDegsConversionFunction')&&(isa(theLattice.MMsToDegsConversionFunction,'function_handle')), 'The %s lattice must have a ''MMsToDegsConversionFunction'' function handle field.', latticeName);


            assert(isfield(theLattice, 'RFpositionsMicrons'), 'The %s lattice must have an ''RFpositionsMicrons'' field.', latticeName);
            assert(isfield(theLattice, 'RFspacingsMicrons'), 'The %s lattice must have an ''RFspacingsMicrons'' field.', latticeName);
            
            % Matrix dimensions must be valid
            [nPos,dimensions] = size(theLattice.RFpositionsMicrons);

            assert(dimensions == 2, ...
                'The ''RFpositionsMicrons'' field of the %s lattice must be an N x 2 matrix. The passed data is %d x %d', latticeName, nPos, dimensions);
            [nPos2,dimensions2] = size(theLattice.RFspacingsMicrons);
            if (nPos2 == 1) 
                theLattice.RFspacingsMicrons = theLattice.RFspacingsMicrons(:);
                [nPos2,dimensions2] = size(theLattice.RFspacingsMicrons);
            end

            assert(dimensions2 == 1, ...
                'The ''RFspacingsMicrons'' field of the %s lattice must be an N x 1 matrix. The passed data is %d x %d', latticeName, nPos2, dimensions2);
            assert(nPos == nPos2, 'The ''RFpositionsMicrons'' field of the %s lattice does not have the same rows (%d) as the ''RFspacingsMicrons'' field (%d)', latticeName,nPos, nPos2); 
        end
           
    end

    % Static methods
    methods (Static)
        [f,v] = facesAndVertices(positions, spacings, shapeOutline);
        theSmoothedSpacings = smoothSpacings(rfSpacings, nearbyRFindices);
        [D,idx] = pdist2(A, B, varargin);
    end

end



