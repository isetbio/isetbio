classdef rgcLayer < handle
    %RGCLAYER - Create an rgcLayer object
    %
    % The rgcLayer class stores general properties of the RGC layer patch
    % and stores the rgcMosaic objects in its mosaic property field.  (This
    % object replaces the innerretina object.)
    %
    %   obj = rgcLayer(inputObj, params);    
    %
    % Usually called internally from rgcLayerCreate ... not sure why. Maybe
    % that will change?  It seems like the input parameters to rgcLayer()
    % should be those of the 
    %
    % An ir object takes as input a bipolar object or an outerSegment object.
    % The ir (inner retina) object stores basic properties about the inner
    % retina such as the position of the simulated retinal patch.
    %
    % See Pillow, Jonathan W., et al. "Spatio-temporal correlations and visual
    % signalling in a complete neuronal population." Nature 454.7207 (2008)
    % and Chichilnisky, E. J., and Rachel S. Kalmar. "Functional asymmetries
    % in ON and OFF ganglion cells of primate retina." The Journal of
    % Neuroscience 22.7 (2002).
    %
    % Properties:
    %
    %  Established by constructor parameters
    %     name:      animal, ir; example: 'macaque ir'
    %     numberTrials: number of trials for spike generation
    %
    %   Inherited from bipolar input
    %     row:       N Stimulus row samples
    %     col:       N Stimulus col samples
    %     size:      Stimulus input spacing (m)
    %     timing:    Stimulus input time step (sec)
    %
    %  Established by the mosaicCreate method
    %     mosaic: cell array of rgc mosaics 
    %
    % Methods: 
    %   set, get, compute, plot
    %
    % Examples:
    %   bpL  = bipolarLayer(coneMosaic);
    %   rgcL = rgcLayer(bpLayer,'name','myRGC');
    %
    %   params.name = 'Macaque inner retina 1';
    %   params.eyeSide = 'left'; params.eyeRadius = 2;
    %   rgcL = rgcLayer(bpLayer, params);
    %
    %  ISETBIO wiki: <a href="matlab:
    %  web('https://github.com/isetbio/isetbio/wiki/Retinal-ganglion-cells','-browser')">RGCS</a>.
    %
    % (c) isetbio team
    %
    % 9/2015 JRG
    % 7/2016 JRG updated
    
    %%
    % Public read/write properties
    properties
    end
    
    % Public, read-only properties.
    properties (SetAccess = public, GetAccess = public)
        
        %NAME Name of this innerRetina
        name;
        
        %NUMBERTRIALS Number of trials when computing
        numberTrials;  
        
        %MOSAIC Cell array containing ganglion cell mosaics
        mosaic;        % The spatial sampling differs for each mosaic
                       
    end
    
    % Protected properties; Methods of the parent class and all of its
    % subclasses can set these.
    properties (SetAccess = protected)
        % Inherited from prior stages
        
        % human, macaque and someday other stuff like mouse
        species;
        
        %SIZE Patch size (m) measured at the cone mosaic (height, width)
        size;        
        
        %TIMESTEP Stimulus temporal sampling (sec) from bipolar
        timeStep;   % This is the same for all mosaics
        
        %EYESIDE Left or right eye
        eyeSide;           
        
        % CENTER of patch (m) with respect to fovea = [0,0];
        center;
        
        % INPUT bipolar layer
        input;
        
    end
    
    % Private properties. Only methods of the parent class can set these
    properties(Access = private)
    end
    
    % Public methods
    methods
        function obj = rgcLayer(bp, varargin)
            % Constructor
            %
            % Taken from the irCreate() code, which was how we called the
            % innerretina object.  We should get rid of this rgcLayerCreate
            % analog and
            %
            %  obj = rgcLayer(bipolarLayer, params)
            %  
            % Required Inputs:
            %  bpLayer:   a bipolar layer object
            %
            % Optional inputs
            %  name:      Name for this instance
            %  species:  
            %  nTrials
            %
            % Outputs:
            %  rgc layer object
            %
            % The coupled-GLM model is described in Pillow, Shlens, Paninski, Sher,
            % Litke, Chichilnisky & Simoncelli, Nature (2008).
            %
            % The LNP and GLM models here are based on
            % <http://pillowlab.princeton.edu/code_GLM.html code by Pillow>
            % under the GNU General Public License.
            %
            % Example:
            %   os  = osCreate('identity');
            %   innerRetina = irCreate(os,'GLM','name','myRGC');
            %   innerRetina = irCreate(os,'name','EJ',...
            %             'eyeSide','left','eyeRadius',12,'eyeAngle',90));
            %
            % See also:  ir.m, rgcMosaic.m
            %
            % JRG 9/2015 Copyright ISETBIO Team
            % JRG 7/2016 updated
            
            % parse input
            p = inputParser;
            
            % Should this by a bipolarLayer??
            p.addRequired('inputObj',@(x)(isa(bp,'bipolarLayer')));
            
            p.addParameter('name','ir1',@ischar);
            p.addParameter('species','macaque',@ischar);
            p.addParameter('nTrials',1,@isscalar);
            
            p.KeepUnmatched = true;
            
            p.parse(bp,varargin{:});
            obj.name         = p.Results.name;
            obj.species      = p.Results.species;
            obj.numberTrials = p.Results.nTrials;
            
            % Should match the cone mosaic patch size and time step
            obj.eyeSide   = bp.input.whichEye; % Maybe not needed?
            obj.size      = bp.size;        % Bipolar patch size
            obj.timeStep  = bp.timeStep;    % Temporal sampling
                        
            % Empty cell array to hold mosaics we create.
            obj.mosaic = cell(1); % Cells are added by mosaicCreate method
            
            % Spatial position on the retina (meters, fovea is 0,0).
            obj.center = bp.center;
            obj.input  = bp;   % Bipolar layer link kept here
            
        end
        
        function obj = mosaicCreate(varargin)
            obj = rgcMosaicCreate(varargin{:});
        end
        
        % set function, see irSet
        function obj = set(obj, varargin)
            obj = irSet(obj, varargin{:});
        end
        
        % get function, see irGet
        function val = get(obj, varargin)
            val = irGet(obj, varargin{:});
        end
        
        %         % IR Compute functions, that loop over the rgc mosaics
        %         function [obj, nTrialsSpikes] = compute(obj, inputObj, varargin)
        %             [obj, nTrialsSpikes] = irCompute(obj,  inputObj, varargin{:});
        %         end
        %
        %         function obj = computeLinearSTSeparable(obj,varargin)
        %             obj = irComputeLinearSTSeparable(obj,varargin{:});
        %         end
        %
        %         function obj = computeSpikes(obj, varargin)
        %             obj = irComputeSpikes(obj,  varargin{:});
        %         end
        %
        %         % plot function, see irPlot
        %         function plot(obj, varargin)
        %             irPlot(obj, varargin{:});
        %         end
        %
        %         % normalize function, see irNormalize
        %         function obj = normalize(obj,varargin)
        %             obj = irNormalize(obj, varargin{:});
        %         end
        
        function val = eccentricity(obj,varargin)
            % Default is units of meters.
            p = inputParser;
            
            % Should check for valid units
            p.addParameter('units','m',@ischar);
            p.parse(varargin{:});
            
            units = p.Results.units;
            val = sqrt(sum(obj.center.^2));
            val = val*ieUnitScaleFactor(units);
            
        end
        
    end
    
    % Methods that must only be implemented in the subclasses. 
    methods (Abstract, Access=public)

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


