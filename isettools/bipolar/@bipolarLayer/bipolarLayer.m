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
        %NUMBERTRIALS Number of trials when computing
        numberTrials;  
        
        %MOSAIC Cell array containing ganglion cell mosaics
        mosaic;        % The spatial sampling differs for each mosaic
                       
    end
    
    % Protected properties; Methods of the parent class and all of its
    % subclasses can set these.
    properties (SetAccess = protected)
        %NAME Name of this innerRetina
        name;       % Note: The computation is specified by the ir subclass
                    % Is the spatial sampling is determined by the bipolar
                    % input?
                    
        %ROW N Stimulus row samples (from bipolar)
        row;        
        
        %COL N Stimulus col samples (from bipolar)
        col;         
        
        %SIZE Patch size (m) measured at the cone mosaic
        size;        
        
        %TIMESTEP Stimulus temporal sampling (sec) from bipolar
        timeStep;   % This is the same for all mosaics
        
        %EYESIDE Left or right eye
        eyeSide;           
        
        %EYERADIUS Position of patch in radius
        eyeRadius;         
        
        %EYEANGLE and angle (degrees)
        eyeAngle;         
        
        %TEMPORALEQUIVECC Temporal equivalent eccentricity (mm)
        temporalEquivEcc; 
        
    end
    
    % Private properties. Only methods of the parent class can set these
    properties(Access = private)
    end
    
    % Public methods
    methods
        function obj = bipolarLayer(bp, varargin)
            % Constructor
            %
            % Taken from the irCreate() code, which was how we called the
            % innerretina object.  We should get rid of this rgcLayerCreate
            % analog and
            %
            %  obj = irCreate(inputObj, params)
            %
            %           params:
            %             'name', name
            %             'eyeSide',{'left','right'},
            %             'eyeRadius',eyeRadius
            %             'eyeAngle',eyeAngle
            %
            % Inputs:
            %  inputObj:  a bipolar object or osDisplayRGBo object
            %  name:      Name for this instance
            %  model:     Computation type for all the rgc mosaics
            %       'linear' - Basic linear filtering as in typical center-surround
            %       'LNP'    - linear-nonlinear-Poisson, see Pillow paper as below; only contains post-spike filter
            %       'GLM'    - coupled generalized linear model, see Pillow paper, includes coupling filters
            %       'Phys'   - pulls a set of parameters measured in physiology by
            %                           the Chichilnisky lab.
            %
            % Outputs:
            %  rgc object: of the specified model type
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

            % We require an inner retina to receive its inputs from a
            % bipolar object.  To skip the bipolar model use a bpIdentity
            % object.
            
            % parse input
            p = inputParser;
            p.addRequired('inputObj',@(x)(isa(bp,'bipolar')||isa(bp{1},'bipolar')));
            
            p.addParameter('name','ir1',@ischar);
            p.addParameter('eyeSide','left',@ischar);
            p.addParameter('eyeRadius',0,@isnumeric);
            p.addParameter('eyeAngle',0,@isnumeric);
            p.addParameter('species','macaque',@ischar);
            p.addParameter('nTrials',1,@isscalar);
            
            p.KeepUnmatched = true;
            
            p.parse(bp,varargin{:});
            
            obj.eyeSide   = p.Results.eyeSide;
            obj.eyeRadius = p.Results.eyeRadius;
            obj.eyeAngle  = p.Results.eyeAngle;
            obj.name      = p.Results.name;
            
            obj.numberTrials = p.Results.nTrials;
            
            if length(bp) > 1                
                obj.size      = bp{1}.get('patch size'); % Bipolar patch
                obj.timeStep  = bp{1}.get('time step');  % Temporal sampling
                
                bpC = bp{1}.get('bipolarResponseCenter');
            else
                obj.size      = bp.get('patch size'); % Bipolar patch
                obj.timeStep  = bp.get('time step');  % Temporal sampling
                
                bpC = bp.get('bipolarResponseCenter');
            end
            
            obj.row = size(bpC,1);  obj.col = size(bpC,2);

            obj.mosaic = cell(1); % Cells are added by mosaicCreate method
            
            % Temporal equivalent eccentricity in deg
            obj.temporalEquivEcc = retinalLocationToTEE(obj.eyeAngle, obj.eyeRadius, obj.eyeSide);
            
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
        
        % IR Compute functions, that loop over the rgc mosaics
        function [obj, nTrialsSpikes] = compute(obj, inputObj, varargin)
            [obj, nTrialsSpikes] = irCompute(obj,  inputObj, varargin{:});
        end
        
        function obj = computeLinearSTSeparable(obj,varargin)
            obj = irComputeLinearSTSeparable(obj,varargin{:});
        end
        
        function obj = computeSpikes(obj, varargin)
            obj = irComputeSpikes(obj,  varargin{:});
        end
        
        % plot function, see irPlot
        function plot(obj, varargin)
            irPlot(obj, varargin{:});
        end
        
        % normalize function, see irNormalize
        function obj = normalize(obj,varargin)
            obj = irNormalize(obj, varargin{:});
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


