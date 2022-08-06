% Subclass of @MosaicConnector, for connecting mRGClattice to a cone lattice
%

classdef coneToMidgetRGCConnector < MosaicConnector

    % Public properties (specific to the coneToMidgetRGCConnector) 
    properties
        
    end
    
    % --- PRIVATE PROPERTIES ----------------------------------------------
    properties (Access = private)              

    end
    % --- END OF PRIVATE PROPERTIES ---------------------------------------
    
    
    % Public methods
    methods
        % Constructor
        function obj = coneToMidgetRGCConnector(...
                sourceLattice, destinationLattice, varargin) 
            p = inputParser;
            p.addParameter('verbosity', 1);
            p.addParameter('coneIndicesToBeConnected', []);
            % Execute the parser
            p.parse(varargin{:});


            % Call the super-class constructor.
            obj = obj@MosaicConnector(...
                sourceLattice, ...
                destinationLattice, ...
                'verbosity', p.Results.verbosity, ...
                'connectableSourceRFindices', p.Results.coneIndicesToBeConnected);
           
        end % Constructor
    end % Public methods

    % Implementations of required -- Public -- Abstract methods defined in the MosaicConnector interface   
    methods
        % mRGCmosaic-specific method to crop the source lattice
        % (depending on the destination lattice)
        cropSourceLattice(obj);

        % mRGCmosaic-specific-specific method to crop the destination lattice
        % (depending on the source lattice)
        cropDestinationLattice(obj);
    end % Implementations of required -- Public -- Abstract methods defined in the MosaicConnector interface
    

    methods (Access = private) 
    end % Private methods 

end


