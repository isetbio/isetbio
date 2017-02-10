classdef rgcLNP < rgcMosaic
% RGCLNP - rgcMosaic using a LNP (linear-nonlinear-Poisson) computational model
% rgcLNP is a subclass of rgcMosaic. It is called when creating a new LNP
% model rgcMosaic for an inner retina object.  Typically we get here from
% the inner retina object with the call:
%
%   ir.mosaicCreate('model','LNP','type','your type goes here')
% 
% The LNP model is detailed in Chichilnisky & Kalmar, J. Neurosci (2002);
% Pillow, Paninski, Uzzell, Simoncelli & Chichilnisky, J. Neurosci (2005);
% and Pillow, Shlens, Paninski, Sher, Litke, Chichilnisky & Simoncelli,
% Nature (2008).
% 
% The computational model implemented here relies on code by
%  the: <a href="matlab:
%  web(http://pillowlab.princeton.edu/code_GLM.html','-browser')">Pillow Lab</a>.
% , which is
% distributed under the GNU General Public License.
% 
% See also: rgcMosaic.m, rgcLinear.m, rgcGLM.m
%
% Example:
% 
%   ir.mosaicCreate('model','LNP','type','on midget'); 
%
% 9/2015 JRG (c) isetbio team
% 7/2016 JRG updated

%% Properties
    % Public, read-only properties.
    properties (SetAccess = public, GetAccess = public)
        
        % Should this be in rgcMosaic or in each of the special types of
        % rgcMosaics? (BW)
        
        % NUMBERTRIALS Number of repeats
        numberTrials = 10;
        
        % Parameter to specify the time bins Pillow uses for coupling and
        % post spike filters (10 ms default)
        dt = 0.1;
    end
           
    % Protected properties.
    properties (SetAccess = private, GetAccess = public)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % rgcLNP is a subclass of rgcMosaic.
        % See the properties in rgcMosaic 
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%              
               
        % GENERATORFUNCTION Pillow promotes the linear input voltage using a nonlinear
        % function that he calls the generator function.  By default this
        % is an exponential.
        generatorFunction;       
        
        % RESPONSEVOLTAGE The nonlinear voltage response after application of the generator
        % function and the spike coupling responses is represented here
        responseVoltage;

        % POSTSPIKEFILTER These hold the parameters used in the computation.
        % This is the response after a spike
        postSpikeFilter;
              
    end
    
    % Private properties. Only methods of the parent class can set these
    properties(Access = private)
    end
    
    % Public methods
    methods
        
        % Constructor
        function obj = rgcLNP(rgc, mosaicType)
            
            % Initialize the parent class            
            obj = obj@rgcMosaic(rgc, mosaicType);

            % Effect of a spike on output voltages
            obj.generatorFunction = @exp;
            
            % Post spike filter
            obj.postSpikeFilter = buildPostSpikeFilter(obj.dt);
            
        end
        
        % set function, see @rgcLNP/mosaicSet.m for details
        function obj = set(obj, varargin)
            mosaicSet(obj, varargin{:});
        end
        
        % get function, see @rgcLNP/mosaicSet.m for details
        function val = get(obj, varargin)
           val = mosaicGet(obj, varargin{:});
        end
      
    end
    
    % Methods that must only be implemented (Abstract in parent class).
    methods (Access=public)               
        initialize(obj);
    end    
    
    % Methods may be called by the subclasses, but are otherwise private 
    methods (Access = protected)
    end
    
    % Methods that are totally private (subclasses cannot call these)
    methods (Access = private)

    end
    
end
