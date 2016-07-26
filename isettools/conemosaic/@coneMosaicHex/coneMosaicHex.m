classdef coneMosaicHex < coneMosaic
    % Create a hexagonal cone mosaic class
    %
    %   cMosaicHex =  coneMosaicHex(resamplingFactor, varargin);
    %
    % The cone mosaic hex defines an array of cones placed on a hexagonal grid
    % which is sampled according to the resamplingFactor. 
    %
    % NPC ISETBIO Team, 2016
    
    properties (SetAccess=private)
        resamplingFactor                % the resamplingFactor used
        coneLocsHexGrid                 % computed coneLocs (hex grid)
        originalResConeLocs             % coneLocs of the originating rect grid
        originalResPattern              % cone pattern of the originaing rect grid
        originalResPatternSampleSize    % original pattern sample size
        originalResFOV                  % original FOV
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
        [activationImage, activationImageLMScone, imageXaxis, imageYaxis] = computeActivationImage(obj, activation);
        
        % Visualize activation images for the hex mosaic (all cones +  LMS submosaics)
        visualizeActivation(obj, activation);
        
        % Print various infos about the cone mosaic
        displayInfo(obj);
        
    end % Public methods
    
    methods (Access = private)
        saveOriginalResState(obj);
        restoreOriginalResState(obj);
    end % Private methods
    
end

