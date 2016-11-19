classdef oiSequence
    % Class for generating a sequence of optical images
    %
    % Usage:
    % (1) theOIsequence = oiSequence(oiBackground, oiModulated, modulationFunction)
    %
    % (2) modulationRegion.radiusInMicrons = 300;
    %     theOIsequence = oiSequence(oiBackground, oiModulated, modulationFunction, 'modulationRegion', modulationRegion, 'oiModulatedReplacesBackground', true);
    %
    % TODO:  Maybe set up some oiGet/Set routines with the syntax
    %        oiSequence.get('oiFixed mumble') and
    %        oiSequence.get('oiModulated mumble')???
    %        Maybe oiSequence.get('frame',val)???
    %
    % See t_oiSequence for example usage.
    %
    %  NPC, ISETBIO TEAM, 2016
    
    properties
        
    end
    
    properties (Dependent)
        length
    end
    
    properties (SetAccess = private)
        % the fixed oi (an oi, the background)
        oiFixed
        
        % the modulated oi (an oi, the pedestal)
        oiModulated
        
        % the oiSequence timebase
        oiTimeAxis
        
        % whether to add the oiModulated to the oiFixed or to blend it with the oiFixed
        composition;
        
        % the modulating function (an array of modulation values, one for
        % each frame)
        modulationFunction;
        
        % the modulating region (a struct describing the region extent to
        % be modulated (for now just a radius)
        modulationRegion;
    end
    
    
    methods  % public methods
        
        % constructor
        function obj = oiSequence(oiFixed, oiModulated, oiTimeAxis, modulationFunction, varargin)
            
            defaultModulationRegion = struct(...
                'radiusInMicrons', nan);
            
            p = inputParser;
            p.addRequired('oiFixed',  @isstruct);
            p.addRequired('oiModulated',  @isstruct);
            p.addRequired('oiTimeAxis', @isnumeric);
            p.addRequired('modulationFunction',  @isnumeric);
            p.addParameter('modulationRegion', defaultModulationRegion, @isstruct);
            p.addParameter('composition', 'add', @ischar);
            p.parse(oiFixed, oiModulated, oiTimeAxis, modulationFunction, varargin{:});
            
            obj.oiFixed = p.Results.oiFixed;
            obj.oiModulated = p.Results.oiModulated;
            obj.oiTimeAxis = p.Results.oiTimeAxis;
            obj.modulationFunction = p.Results.modulationFunction;
            obj.modulationRegion = p.Results.modulationRegion;
            obj.composition = p.Results.composition;
            
            % The oiTimeAxis must be the same length as the number of
            % values in the modulationFunction.  If it only has 1 value, we
            % are going to assume that the value is delta T and we will
            % create the whole vector.  If it has a vector, but that vector
            % is not the same length as the modulationFunction, we throw an
            % error.
            if length(obj.oiTimeAxis) == 1
                % First moment in time is 0. Increments by the set value.
                obj.oiTimeAxis = obj.oiTimeAxis*(0:(length(modulationFunction))-1);
            elseif length(obj.oiTimeAxis) ~= length(obj.modulationFunction)
                error('Time axis does not match modulation function');
            end
            
            % Set a validation function above, don't do this.
            if (~strcmp(obj.composition, 'add')) && (~strcmp(obj.composition, 'blend'))
                error('''composition'' must be set to either ''blend'' or ''add''.');
            end
            
            % Make sure that oiFixed and oiModulated have identical shape
            oiFixedSpatialSupport      = round(oiGet(obj.oiFixed, 'spatial support','microns'), 7);
            oiModulatedSpatialSupport  = round(oiGet(obj.oiModulated, 'spatial support','microns'), 7);
            
            if (any(size(oiFixedSpatialSupport) ~= size(oiModulatedSpatialSupport)))
                error('Mismatch between spatial dimensions of oiFixed, oiModulated');
            end
            if (any(oiFixedSpatialSupport(:) ~= oiModulatedSpatialSupport(:)))
                error('Mismatch between spatial support of oiFixed, oiModulated');
            end
        end
        
        %% Define methods in the @oiSequence directory
        
        % Method for on-the-fly computation of the oi at desired index
        oiFrame = frameAtIndex(obj, index);
        
        % Visualize the sequence
        visualize(obj,varargin);
        
        % Use the eye movements?  We have eye movements in coneMosaic, not
        % oi.  So, I think we need to wait to see the eye movements until
        % we are in the cMosaic case.  If we want them here, then we need
        % to take a whole em object. (BW).
        % visualizeWithEyeMovementSequence(obj, emTimeAxis);
        
        %% Local get methods
        
        % Return the length of the oiSequence
        function val = get.length(obj)
            val = numel(obj.modulationFunction);
        end
        
        % Return the modulationFunction used
        function val = get.modulationFunction(obj)
            val = obj.modulationFunction;
        end
        
        % Return the oiTimeAxis used
        function val = get.oiTimeAxis(obj)
            val = obj.oiTimeAxis;
        end
        
        % Return the composition type used
        function val = get.composition(obj)
            val = obj.composition;
        end
        
        
    end
    
end

