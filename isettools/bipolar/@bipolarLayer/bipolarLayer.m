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
    % Optional Key/Value Pairs:
    %    name           - the name of the bipolarLayer object.
    %                     (Default 'bipolarLayer')
    %    nTrials        - the number of trials. (Default 1)
    %
    % Outputs:
    %    bipolarLayer   - A bipolarLayer object
    %

    % History:
    % BW (c) isetbio team
    %
    %    10/18/17  jnm  Comments & formatting

    %% Public Properties
    % Public read/write properties
    properties
    end
    %%%
    % Public, read-only properties.
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
        function obj = bipolarLayer(cMosaic, varargin)
            %% Create a Bipolar Layer object
            %
            % Taken from the irCreate() code, which was how we called the
            % innerretina object.  We should get rid of this rgcLayerCreate
            % analog and
            %
            %  obj = bipolarLayer(coneMosaic, params)
            %
            %  Required:
            %    coneMosaic object
            %
            %  Optional
            %   'name'    -   string
            %   'nTrials' -   1
            %
            % Outputs:
            %  bipolarLayer 
            %
            % BW, ISETBIO Team, 2017
            %%%
            % parse input
            p = inputParser;
            p.addRequired('cMosaic',@(x)(isa(cMosaic,'coneMosaic')));
            
            p.addParameter('name','bipolarLayer',@ischar);
            p.addParameter('nTrials',1,@isscalar);
            
            p.KeepUnmatched = true;
            
            p.parse(cMosaic,varargin{:});
            
            obj.name         = p.Results.name;
            obj.nTrials = p.Results.nTrials;

            % We may keep the cone mosaic around and then get rid of the
            % obj.center, .species, .timeStep and .size.  No reason to have
            % both other than confusion.  We can write little methods to
            % get these from the input.
            obj.input  = cMosaic;

            % The patch size should apply to all of the mosaics.  It is
            % determined by the cone mosaic patch, as is the time step.
            obj.size = cMosaic.size;
            %%%
            % The time sample is the integration time of the cones.  Both
            % the absorptions and the current are sampled at this rate
            obj.timeStep  = cMosaic.integrationTime;  

            % Center of the patch with respect to distance (meters) on the
            % retina.  Fovea is (0,0).
            obj.center =  cMosaic.center;
        end

        % Return Window handle
        function hdl = window(obj,varargin)
            % Show the bipolar layer window
            hdl = bipolarLayerWindow(obj);
        end
    end
    
    %% Methods that must only be implemented in the subclasses. 
    methods (Access=public)
        %% Public Subclass Methods
        % Could define the methods here, or just let them be implicitly
        % defined in the functions themselves
        hdl = plot(obj,pType,varargin);
        str = describe(obj);
    end
    %% Protected Subclass Methods
    % Methods may be called by the subclasses, but are otherwise private
    methods (Access = protected)
        % 1D spatial convolution of array of center & surround responses
        spConvolve(obj);
        % 1D temporal convolution of array of center & surround responses
        timeConvolve(obj);
    end
    %% Private Subclass Methods
    % Methods that are totally private (subclasses cannot call these)
    methods (Access = private)
    end
end
