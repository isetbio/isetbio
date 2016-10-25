classdef oiSequence
    % Class for generating a sequence of optical images
    %
    % Usage:
    % (1) theOIsequence = oiSequence(oiBackground, oiModulated, modulationFunction)
    %
    % (2) modulationRegion.radiusInMicrons = 300;
    %     theOIsequence = oiSequence(oiBackground, oiModulated, modulationFunction, 'modulationRegion', modulationRegion, 'oiModulatedReplacesBackground', true);
    %
    % See t_oiSequence for example usage.
    %
    %  NPC, ISETBIO TEAM, 2016
    
    properties 
        
    end
    
    properties (Dependent)
        length
    end
    
    properties (Access = private)   
        % the fixed oi (an oi, the background)
        oiFixed
        fixedPhotons
        
        % the modulated oi (an oi, the pedestal)
        oiModulated
        modulatedPhotons
        
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
        function obj = oiSequence(oiFixed, oiModulated, modulationFunction, varargin)
            
            defaultModulationRegion = struct(...
                'radiusInMicrons', nan);
            
            p = inputParser;
            p.addRequired('oiFixed',  @isstruct);
            p.addRequired('oiModulated',  @isstruct);
            p.addRequired('modulationFunction',  @isnumeric);
            p.addParameter('modulationRegion', defaultModulationRegion, @isstruct);
            p.addParameter('composition', 'add', @ischar);
            p.parse(oiFixed, oiModulated, modulationFunction, varargin{:});
            
            obj.oiFixed = p.Results.oiFixed;
            obj.oiModulated = p.Results.oiModulated;
            obj.modulationFunction = p.Results.modulationFunction;
            obj.modulationRegion = p.Results.modulationRegion;
            obj.composition = p.Results.composition;
            
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
            % Extract the photons
            obj.fixedPhotons = oiGet(obj.oiFixed, 'photons');
            obj.modulatedPhotons = oiGet(obj.oiModulated, 'photons');
        end
            
        function val = get.length(obj)
            val = numel(obj.modulationFunction);
        end
        
        function oiFrame = frameAtIndex(obj, frameIndex)

            if (strcmp(obj.composition, 'add'))
                retinalPhotons = obj.fixedPhotons + obj.modulationFunction(frameIndex)*obj.modulatedPhotons;
            else
                retinalPhotons = obj.fixedPhotons*(1-obj.modulationFunction(frameIndex)) + obj.modulationFunction(frameIndex)*obj.modulatedPhotons;
            end
            
            if (~isnan(obj.modulationRegion.radiusInMicrons))
                % modulate a subregion only
                pos = oiGet(obj.oiModulated, 'spatial support', 'microns');
                ecc = sqrt(sum(pos.^2, 3));
                mask = ecc < obj.modulationRegion.radiusInMicrons;
                
                for k = 1:size(retinalPhotons,3)
                    fullFrame = retinalPhotons(:,:, k);
                    background = obj.fixedPhotons;
                    background = background(:,:, k);
                    fullFrame(mask == 0) = background(mask == 0);
                    retinalPhotons(:,:, k) = fullFrame;
                end
            end
            
            % Return the composite oi
            oiFrame = obj.oiFixed;
            oiFrame = oiSet(oiFrame, 'photons', retinalPhotons);
        end
    end
    
end

