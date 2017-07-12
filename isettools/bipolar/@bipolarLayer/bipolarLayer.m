classdef bipolarLayer < handle
    % BIPOLARLAYER - Create an bipolarLayer object
    %
    % The bipolarLayer class stores general properties of the bipolar layer
    % patch and stores the bipolarMosaic objects in its mosaic property
    % field. This object has a role that matches the rgcLayer object.
    %
    %   obj = bipolarLayer([conemosaicObject], params);    
    %
    % Properties:
    %
    % BW (c) isetbio team
    
    %%
    % Public read/write properties
    properties
    end
    
    % Public, read-only properties.
    properties (SetAccess = public, GetAccess = public)
        %NAME Name of this bipolar layer object
        % Could be used in window 
        name;       
        
        %NUMBERTRIALS Number of trials when computing
        numberTrials;  
        
        %MOSAIC Cell array containing bipolar cell mosaics
        mosaic;
        
        % When we have a window, we use this figureHandle for refresh
        figureHandle;
                       
    end
    
    % Protected properties; Methods of the parent class and all of its
    % subclasses can set these.
    properties (SetAccess = protected)
        % These are protected because they are determined from the cone
        % mosaic that provides the input and thus should not change
        % A few parameters stored here for convenience, but they can be
        % derived from input or input to input or ...
        
        %TIMESTEP Stimulus temporal sampling (sec) from bipolar
        timeStep;   % This is the same for all mosaics 
        
        %CENTER position of the patch with respect to fovea (0,0)
        center;
        
        %SIZE Patch size (m) measured at the cone mosaic
        size;  
        
        % INPUT  - Cone mosaic input
        input;
        
    end
    
    % Private properties. Only methods of the parent class can set these
    properties(Access = private)
    end
    
    % Public methods
    methods
        function obj = bipolarLayer(cMosaic, varargin)
            % Constructor
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

            % parse input
            p = inputParser;
            p.addRequired('cMosaic',@(x)(isa(cMosaic,'coneMosaic')));
            
            p.addParameter('name','bipolarLayer',@ischar);
            p.addParameter('nTrials',1,@isscalar);
            
            p.KeepUnmatched = true;
            
            p.parse(cMosaic,varargin{:});
            
            obj.name         = p.Results.name;
            obj.numberTrials = p.Results.nTrials;
            
            % Create an empty cell array of bipolar mosaics
            obj.mosaic = [];
            
            % We may keep the cone mosaic around and then get rid of the
            % obj.center, .species, .timeStep and .size.  No reason to have
            % both other than confusion.  We can write little methods to
            % get these from the input.
            obj.input  = cMosaic;
            
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
        
        % Show the bipolar layer window
        function hdl = window(obj,varargin)
            hdl = bipolarLayerWindow(obj);
        end

    end
    
    % Methods that must only be implemented in the subclasses. 
    methods (Access=public)
        % Could define the methods here, or just let them be implicitly
        % defined in the functions themselves
        hdl = plot(obj,pType,varargin);
        str = describe(obj);
    end
    
    % Methods may be called by the subclasses, but are otherwise private
    methods (Access = protected)
        spConvolve(obj);
        timeConvolve(obj);
    end
    
    % Methods that are totally private (subclasses cannot call these)
    methods (Access = private)
    end
    
end


