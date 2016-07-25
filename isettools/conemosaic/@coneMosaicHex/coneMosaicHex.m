classdef coneMosaicHex < coneMosaic
    % Create a hexagonal cone mosaic class
    %
    %   cMosaicHex =  coneMosaicHex('cone', cone, 'os', os);
    %
    % The cone mosaic defines an array of cones.  The individual cones have
    % absorption properties defined by cMosaic.cone.  The computation from
    % absorptions to photocurrent is defined by cMosaic.os
    %
    % NPC ISETBIO Team, 2016
    
    properties (SetAccess=private)
        upSampleFactor
        coneLocsHexGrid                 % computed coneLocs (hex grid)
        originalResConeLocs             % original res coneLocs (rect grid)
        originalResPattern
        originalResPatternSampleSize
        originalResFOV
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
            p.addRequired('upSampleFactor', @isnumeric);
            p.parse(upSampleFactor);
            obj.upSampleFactor = p.Results.upSampleFactor;
            
            % Generate sampled hex grid
            obj.resampleGrid(obj.upSampleFactor);
        end
        
        setSizeToFOVForHexMosaic(obj,fov);
        resampleGrid(obj, upSampleFactor);
        visualizeGrid(obj, varargin);
        visualizeActivation(obj, activation);
        [activationImage, imageXaxis, imageYaxis] = computeActivationImage(obj, activation);
        displayInfo(obj);
        
    end % Public methods
    
    methods (Access = private)
        computeHexMosaic(obj);
        saveOriginalResState(obj);
        restoreOriginalResState(obj);
    end % Private methods
    
end

