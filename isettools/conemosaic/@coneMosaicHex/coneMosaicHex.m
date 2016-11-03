classdef coneMosaicHex < coneMosaic
    % Create a hexagonal cone mosaic class
    %
    %   cMosaicHex =  coneMosaicHex(resamplingFactor, varyingDensity, varargin);
    %
    % The cone mosaic HEX is a subclass of coneMosaic. It differs because
    % the array of cones is placed on a hexagonal, rather than rectangular,
    % grid. 
    %
    % The hex mosaic is sampled according to the resamplingFactor. The cone
    % density can be spatially-varying if varyingDensity is set to true.
    %
    % Name-Value parameter pairs used with coneMosaic also can be used with
    % coneMosaicHex, e.g., 
    %
    % Example:
    %      resamplingFactor = 8;
    %      varyingDensity = false;
    % cMosaicHex = coneMosaicHex(resamplingFactor, varyingDensity, ...
    %                          'name', 'the hex mosaic', ...
    %                          'size', [48 32], ...
    %                     'noiseFlag', 0,  ...
    %                'spatialDensity', [0 0.6 0.3 0.1] ...
    %       );
    %   cMosaicHex.window;
    %
    % See also: coneMosaic.
    %
    % NPC ISETBIO Team, 2016
    
    properties (SetAccess=private)
        lambdaMin                               % min cone separation in the mosaic
        lambdaMid                               % the cone separation at the middle of the mosaic
        varyingDensity                          % whether to have an eccentricity-based spatially-varying density (boolean)
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
        function obj = coneMosaicHex(upSampleFactor, varyingDensity, varargin)
            % Initialize the hex cone mosaic class
            %   cMosaic =  coneMosaicHex(upSampleFactor, varyingDensity, ['cone',cone,'os','os]);
            
            % Call the super-class constructor.
            obj = obj@coneMosaic(varargin{:});
            
            % Get a copy of the original coneLocs
            obj.saveOriginalResState();
            
            % parse input
            p = inputParser;
            p.addRequired('resamplingFactor', @isnumeric);
            p.addRequired('varyingDensity', @islogical);
            
            p.parse(upSampleFactor, varyingDensity);
            obj.resamplingFactor = p.Results.resamplingFactor;
            obj.varyingDensity = p.Results.varyingDensity;
            
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
        hFig = visualizeActivationMaps(obj, activation, varargin);
        
        % Print various infos about the cone mosaic
        displayInfo(obj);
        
    end % Public methods
    
    methods (Access = private)
        % Private methods
        %
        % These are templated here. The actual code is in the
        % @coneMosaicHex directory.
        saveOriginalResState(obj);
        restoreOriginalResState(obj);
    end 
    
end

