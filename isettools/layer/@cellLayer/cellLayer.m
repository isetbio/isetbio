classdef cellLayer < handle
    % CELLLAYER - Define a cellLayer object
    %
    % The cellLayer class is the base for bipolarLayer and rgcLayer and
    % probably others in the future. It stores general properties and
    % methods of layers
    %
    %   obj = cellLayer(varargin);
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
        nTrials;
        
        % When we have a layer window, this is the figure handle
        % gdata = guidata(fig); Should work.
        fig;
        
        % NOTES - cell array to store notes
        notes = {};
    end
    
    % Protected properties; Methods of the parent class and all of its
    % subclasses can set these.
    properties (SetAccess = protected)
        % These are protected because they are determined from the cone
        % mosaic that provides the input and thus should not change
        % A few parameters stored here for convenience, but they can be
        % derived from input or input to input or ...
        
        % Cell array of inputs to the mosaics in this layer
        input;
        
        %TIMESTEP Stimulus temporal sampling (sec) from original layer
        % This is the same for all mosaics and thus shared here.
        timeStep;
        
        %CENTER position of the input layers with respect to fovea (0,0)
        center;
        
        %SIZE Patch size (m) measured at the cone mosaic, inherited up the
        %layer chain
        size;
        
    end
    
    % Private properties. Only methods of the parent class can set these
    properties(Access = private)
    end
    
    % Public methods
    methods
        function obj = cellLayer(varargin)
            % Constructor
            %
            % BW, ISETBIO Team, 2017
            
            % parse input
            p = inputParser;
            p.KeepUnmatched = true;
            
            p.addParameter('name','cellLayer',@ischar);
            p.addParameter('nTrials',1,@isscalar);
            
            p.parse(varargin{:});
            
            % We may over-write these with the subclass of layer
            obj.name    = p.Results.name;
            obj.nTrials = p.Results.nTrials;
            
        end
        
        
    end
end
    
