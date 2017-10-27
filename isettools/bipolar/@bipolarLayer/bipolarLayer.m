classdef bipolarLayer < cellLayer
    % Create an bipolarLayer object
    %
    % Syntax:
    %   obj = bipolarLayer([conemosaicObject], params);    
    %
    % Description:
    %    Create an bipolarLayer object
    %
    %    This is a subclass of the general class 'cellLayer'.  It inherits
    %    the properties of that class, including slots for an input, fig,
    %    center, size, ...
    %
    %    The bipolarLayer class stores general properties of the bipolar
    %    layer patch and stores the bipolarMosaic objects in its mosaic
    %    property field. This object has a role that matches the rgcLayer
    %    object.
    %
    % Inputs:
    %    coneMosaic - cone mosaic object including photocurrent response
    %
    % Optional key/value pairs:
    %    'name'           - the name of the bipolarLayer object.
    %                     (Default 'bipolarLayer')
    %    'nTrials'        - the number of trials. (Default 1)
    %
    % Outputs:
    %    obj - Returned bipolarLayer object
    %
    %
    % Notes:
    % * [NOTE: DHB - I started a header comment for the window method, but
    %    this needs to be filled out.]
    %
    % Known bugs/limitations:
    % * [NOTE: DHB - I can find no evidence that methods spConvolve and
    %    timeConvolve, both declared here, are implemented anywhere.
    %    Delete the declarations?]
    %
    % * [NOTE: DHB - Why does public window method show up as protected 
    %    when we use "doc bipolarLayer"?  An ineritance issue with this
    %    method?]
    %
    
    % History:
    %    xx/xx/xx  baw (c) isetbio team
    %    10/18/17  jnm Comments & formatting

    %% Public Properties
    % Public read/write properties
    properties
    end
   
    %% Public, read-only properties.
    properties (SetAccess = public, GetAccess = public)
        % MOSAIC Cell array holding the bipolar cell mosaics
        mosaic;
    end
    
    %% Protected Properties
    % Protected properties; Methods of the parent class and all of its
    % subclasses can set these.
    properties (SetAccess = protected)       
    end
    
    %% Private Properties
    % Private properties. Only methods of the parent class can set these
    properties(Access = private)
    end
    
    %% Public methods
    methods
        % Constructor
        function obj = bipolarLayer(cMosaic, varargin)

            % parse input
            p = inputParser;
            p.addRequired('cMosaic',@(x)(isa(cMosaic,'coneMosaic'))); 
            p.addParameter('name','bipolarLayer',@ischar);
            p.addParameter('nTrials',1,@isscalar);   
            p.KeepUnmatched = true; 
            p.parse(cMosaic,varargin{:});
            
            obj.name = p.Results.name;
            obj.nTrials = p.Results.nTrials;

            % We may keep the cone mosaic around and then get rid of the
            % obj.center, .species, .timeStep and .size.  No reason to have
            % both other than confusion.  We can write little methods to
            % get these from the input.
            obj.input  = cMosaic;

            %%%
            % The patch size should apply to all of the mosaics.  It is
            % determined by the cone mosaic patch, as is the time step.
            obj.size = cMosaic.size;

            % The time sample is the integration time of the cones.  Both
            % the absorptions and the current are sampled at this rate
            obj.timeStep  = cMosaic.integrationTime;  

            % Center of the patch with respect to distance (meters) on the
            % retina.  Fovea is (0,0).
            obj.center =  cMosaic.center;
        end

        function hdl = window(obj,varargin)
            % Return handle to bipolar layer window
            %
            % Syntax:
            %    hdl = obj.bipolarLayerWindow; 
            hdl = bipolarLayerWindow(obj);
            
        end
    end
    
    %% Public Subclass Methods
    %
    % These are methods implemented for this subclass, with
    % their calling form defined here.
    methods (Access=public)
        % % See bipolarLayer/plot.m 
        hdl = plot(obj,pType,varargin);
        
        % % See bipolarLayer/describe.m 
        str = describe(obj);
    end
    
    %% Protected Subclass Methods
    % 
    % These are methods may be called by subclasses, but are otherwise private
    methods (Access = protected)
        % 1D spatial convolution of array of center & surround responses
        spConvolve(obj);
        
        % 1D temporal convolution of array of center & surround responses
        timeConvolve(obj);
        
    end
    
    %% Private Subclass Methods
    %
    % These are methods that are totally private to this class (subclasses cannot call these)
    methods (Access = private)
    end
end
