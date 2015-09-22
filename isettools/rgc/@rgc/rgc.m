classdef rgc < handle
% The @rgc parent class is a component of isetbio that is used to compute
% the responses of retinal gnanglion cells from the isetbio @outerSegment 
% object. Subclasses of @rgc include @rgcLNP and @rgcGLM and are used to
% instantiate specific models of RGC processing.
%     
% (c) isetbnio
% 
% 9/2015 JRG
% 
% 

    % Public read/write properties
    properties
    end
    
    % Public, read-only properties.
    properties (SetAccess = private, GetAccess = public)  

    end
    
    % Protected properties; Methods of the parent class and all of its
    % subclasses can set these.
    properties (SetAccess = protected)
        noiseFlag;
    end
    
    % Private properties. Only methods of the parent class can set these
    properties(Access = private)
        
    end
    
    % Public methods
    methods
        
        function obj = rgc()
            obj.initialize();
        end
        
        % set function, see
        function obj = set(obj, param, val, varargin)
            rgcSet(obj, param, val, varargin);
        end
        
        % get function, see 
        function val = get(obj, param, varargin)
           val = rgcGet(obj, param, varargin);
        end
        
    end
    
    % Methods that must only be implemented in the subclasses.
    % If a subclass does not implement each and every of these methods
    % it cannot instantiate objects.
    
    methods (Abstract, Access=public)
        % see 
        compute(obj, sensor, param, varargin);
        % see 
        plot(obj, sensor);
    end
    
    % Methods may be called by the subclasses, but are otherwise private 
    methods (Access = protected)
        spConvolve1D
        spConvolve
        tempConvolve
        fullConvolve
       
    end
    
    % Methods that are totally private (subclasses cannot call these)
    methods (Access = private)
        initialize(obj);
    end
    
end

    
    