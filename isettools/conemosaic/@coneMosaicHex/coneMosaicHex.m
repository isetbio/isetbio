classdef coneMosaicHex < coneMosaic
    % Create a hexagonal cone mosaic class
    %
    %   cMosaicHex =  coneMosaicHex(resamplingFactor, varargin);
    %
    % The cone mosaic hex defines an array of cones placed on a hexagonal grid
    % which is sampled according to the resamplingFactor. All key-value
    % parameter pairs used with coneMosaic can be used with coneMosaicHex.
    % e.g.:
    %   cMosaicHex = coneMosaicHex(resamplingFactor, ...
    %                     'name', 'the hex mosaic', ...
    %                      'size', [48 32], ...
    %                   'pattern', LMSpattern, ...
    %                 'noiseFlag', 0,  ...
    %            'spatialDensity', [0 0.6 0.3 0.1] ...
    %   );
    % See also coneMosaic.
    %
    % NPC ISETBIO Team, 2016
    
    properties (SetAccess=private)
        resamplingFactor                        % resamplingFactor
        coneLocsHexGrid                         % computed coneLocs (hex grid)
        coneLocsOriginatingRectGrid             % coneLocs of the originating rect grid
        patternOriginatingRectGrid              % cone pattern of the originating rect grid
        patternSampleSizeOriginatingRectGrid    % pattern sample size of the originating rect grid
        fovOriginatingRectGrid                  % FOV of the originating rect grid
    end
    
    % Public methods
    methods
        
        % Constructor
        function obj = coneMosaicHex(upSampleFactor, varargin)
            % Initialize the hex cone mosaic class
            %   cMosaic =  coneMosaicHex(upSampleFactor, ['cone',cone,'os','os]);
            
            % Call the super-class constructor.
            obj = obj@coneMosaic(varargin{:});
            
            % Get a copy of the original coneLocs
            obj.saveOriginalResState();
            
            % parse input
            p = inputParser;
            p.addRequired('resamplingFactor', @isnumeric);
            p.parse(upSampleFactor);
            obj.resamplingFactor = p.Results.resamplingFactor;
            
            % Generate sampled hex grid
            obj.resampleGrid(obj.resamplingFactor);
        end
        
        % Change the FOV of the mosaic
        setSizeToFOVForHexMosaic(obj,fov);
        
        % Sample the original rectangular mosaic using a hex grid sampled at the passed resamplingFactor
        resampleGrid(obj, resamplingFactor);
        
        % Visualize different aspects of the hex grid
        visualizeGrid(obj, varargin);
        
        % Compute activation images for the hex mosaic (all cones +  LMS submosaics)
        [activationImage, activationImageLMScone, imageXaxis, imageYaxis] = computeActivationDensityMap(obj, activation);
        
        % Visualize activation maps images for the hex mosaic (all cones +  LMS submosaics)
        visualizeActivationMaps(obj, activation, varargin);
        
        % Print various infos about the cone mosaic
        displayInfo(obj);
        
    end % Public methods
    
    methods (Access = private)
        saveOriginalResState(obj);
        restoreOriginalResState(obj);
    end % Private methods
    
end

