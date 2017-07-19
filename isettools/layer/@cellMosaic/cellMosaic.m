classdef cellMosaic < handle
%CELLMOSAIC - Define the cell mosaic class used by bipolarMosaic, rgcMosaic
%             and possibly others
%
%
%     cellType;                        % there are five different cell types
%     cellLocation;                    % location of bipolar RF center
%     patchSize;                       % size of retinal patch from sensor
%     timeStep;                        % time step of simulation from sensor
%     filterType;                      % bipolar temporal filter type
% 
%  ISETBIO wiki: <a href="matlab:
%  web('https://github.com/isetbio/isetbio/wiki/bipolar','-browser')">bipolar</a>.
%   
% 5/2016 JRG (c) isetbio team

%% Define object
% Public, read-only properties.
properties (SetAccess = private, GetAccess = public)
end

% Protected properties.
properties (SetAccess = protected, GetAccess = public)
    % We need an amplitude for the center and surround
    
    % CELLTYPE (e.g., midget, parasol, onsbc, ....)s
    cellType;                        
    
    %CELLLOCATION location of bipolar RF centers w.r.t. the input samples
    cellLocation;
    
    % PATCHSIZE  cone mosaic referred specification of patch size (meters)
    patchSize;                       
    
    % TIMESTEP time step of simulation inherited from cone mosaic
    timeStep;       
    
    % FILTERTYPE temporal filter type
    filterType; 
    
    % SRFCENTER spatial RF of the center on the input samples
    % Represented w.r.t the input sampling grid
    sRFcenter = [];                       
    
    % SRFSURROUND spatial RF of the surround on the input samples
    % Represented w.r.t the input sampling grid
    sRFsurround = [];                   
    
    % input mosaic - could be cell array of mosaics, I suppose
    input;
    
    % parent - the layer object containing this and other mosaics
    parent;

end

% Private properties. Only methods of the parent class can set these
properties(Access = private)
end

properties (Access = public)
    %FIGUREHANDLE When we open the figure for the mosaic, we store the handle here
    fig;
    
    % NOTES - A cell array to keep notes, parameters, unaccounted for
    % things.  Maybe it should be a struct.
    notes = {};
end

% Public methods
methods
    
    % Constructor
    function obj = cellMosaic(varargin)     
        % Initialize the mosaic object.
        % Details are specified by the subclass, such as bipolarMosaic or
        % rgcMosaic, or all the way to rgcGLM
        %
        % BW ISETBIO Team, 2017
        
%         p = inputParser;
%         p.addParameter('input',[]);
%         p.addParameter('parent',[], @(x)(isequal(class(x),'cellLayer')));
%         p.addParameter('filterType',  1, @isnumeric);
%         
%         p.parse(varargin{:});  
%         
%         % The layer object that this is part of.
%         obj.input     = p.Results.input;
%         obj.parent    = p.Results.parent;
% 
%         % This might be a mistake, but we store the size and time step of
%         % the cone mosaic in this mosaic for easy accessibility.  The
%         % cMosaic itself is stored in the layer object that contains this
%         % mosaic.  Maybe these should just be private variables?
%         obj.patchSize = input.size; 
%         obj.timeStep  = input.integrationTime;
                
    end
    
end

% Methods that must only be implemented (Abstract in parent class).

% Methods may be called by the subclasses, but are otherwise private
methods (Access = protected)
end

% Methods that are totally private (subclasses cannot call these)
methods (Access = private)
end

end