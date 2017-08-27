classdef coneMosaicHex < coneMosaic
    %CONEMOSAICHEX Create a hexagonal cone mosaic class
    %
    %   cMosaicHex =  CONEMOSAICHEX(resamplingFactor,  varargin);
    %
    % The cone mosaic HEX is a subclass of coneMosaic. It differs because
    % the array of cones is placed on a hexagonal, rather than rectangular,
    % grid. 
    %
    % The hex mosaic is sampled according to the resamplingFactor. The cone
    % density can be spatially-varying if eccBasedConeDensity is set to true.
    %
    % The customLambda argument is empty to obtain default performance, but
    % may be set to set the spacing for regularly spaced hexagonal mosaics.
    % The units of customLambda are microns
    %
    % Name-Value parameter pairs used with coneMosaic also can be used with
    % coneMosaicHex 
    %
    % Example:
    %      resamplingFactor = 8;
    %      eccBasedConeDensity = false;
    %      customLambda = 3.0;              %If not passed (or set to []) @coneMosaiHex chooses
    %                                       the cone spacing based on the eccentricity 
    %                                       of the mosaic as determined by the 
    %                                       coneSizeReadData function.
    %                                       If set to a value (specified in microns), 
    %                                       cone spacing is set to that value. Note that
    %                                       if the 'eccBasedConeDensity' param is  set to true, 
    %                                       the 'customLambda' param is ignored.
    %     customInnerSegmentDiameter = 2.5; %If not passed (or set to []) @coneMosaiHex chooses
    %                                        the default pigment.pdWidth, pigment.pdHeight values from
    %                                        its superclass. If it set to a value (specified in microns)
    %                                        @coneMosaicHex sets pigment.pdWidth and pigment.pdheight to
    %                                        sizeForSquareApertureFromDiameterForCircularAperture(customInnerSegmentDiameter)
    %
    % Example usage
    % cMosaicHex = coneMosaicHex(resamplingFactor, ...
    %   'name', 'the hex mosaic', ...
    %   'fovDegs', 0.35, ...
    %   'eccBasedConeDensity', true, ...   
    %   'noiseFlag', 0,  ...
    %   'spatialDensity', [0 0.6 0.3 0.1]);
	%
    %   cMosaicHex.window;
    %
    % See also: CONEMOSAIC, t_coneMosaicHex, t_coneMosaicHexReg
    
    % NPC ISETBIO Team, 2016
    
    properties (SetAccess=private)
        lambdaMin                               % min cone separation in the mosaic
        lambdaMid                               % the cone separation at the middle of the mosaic
        customLambda                            % user-supplied lambda (cone spacing) for regularly spaced mosaics (in microns)
        customInnerSegmentDiameter              % user-supplied inner segment diameter (for a circular aperture, in microns)
        eccBasedConeDensity                     % whether to have an eccentricity-based spatially-varying density (boolean)
        resamplingFactor                        % resamplingFactor
        marginF                                 % factor determining how much more (if >1) or less (<1) mosaic to compute
        sConeMinDistanceFactor                  % min distance between neighboring S-cones (to make the S-cone lattice semi-regular) = f * local cone separation 
        sConeFreeRadiusMicrons                  % radius of S-cone free retina, default: 45 microns, which is 0.15, so S-cone free region = 0.3 degs diameter
        coneLocsHexGrid                         % floating point coneLocs on the hex grid. This is sampled according to the resamplingFactor.
        coneLocsOriginatingRectGrid             % coneLocs of the originating rect grid
        patternOriginatingRectGrid              % cone pattern of the originating rect grid
        patternSampleSizeOriginatingRectGrid    % pattern sample size of the originating rect grid
        fovOriginatingRectGrid                  % FOV of the originating rect grid
        rotationDegs                            % rotation in degrees
        saveLatticeAdjustmentProgression        % flag indicating whether to save the iterative lattice adjustment steps
        latticeAdjustmentSteps                  % 3D array with coneLocsHexGrid at each step of the lattice adjustment
        initialLattice                          % coneLocsHexGrid at iteration 0 (perfect hex grid)
        latticeAdjustmentPositionalToleranceF   % tolerance for whether to decide that there is no move movement
    	latticeAdjustmentDelaunayToleranceF     % tolerance for deciding whether to trigger another Delaunay triangularization
    end
    
    % Public methods
    methods
        
        % Constructor
        function obj = coneMosaicHex(upSampleFactor, varargin)
            % Initialize the hex cone mosaic class
            %   cMosaic =  coneMosaicHex(upSampleFactor, ['varyingDensity', true, 'customLambda', 3, 'customInnerSegmentDiameter'', 2.5, 'cone',cone,'os','os]);
            
            % Params that we want to consume (not pass to our super-class @coneMosaic)
            paramsForConeMosaicHex = {...
                'fovDegs', ...
                'eccBasedConeDensity', ...
                'sConeMinDistanceFactor', ...
                'sConeFreeRadiusMicrons', ...
                'customLambda', ...
                'customInnerSegmentDiameter', ...
                'rotationDegs', ...
                'saveLatticeAdjustmentProgression',...
                'latticeAdjustmentPositionalToleranceF', ...
                'latticeAdjustmentDelaunayToleranceF' ...
                'marginF' ...
                };
            
            % Call the super-class constructor.
            vararginForConeMosaic = {};
            vararginForConeHexMosaic = {};
            for k = 1:2:numel(varargin)
                if (ismember(varargin{k}, paramsForConeMosaicHex))
                    vararginForConeHexMosaic{numel(vararginForConeHexMosaic)+1} = varargin{k};
                    vararginForConeHexMosaic{numel(vararginForConeHexMosaic)+1} = varargin{k+1};
                else
                    if (strcmp(varargin{k}, 'center'))
                        error('Currently coneMosaicHex only supports mosaics centered at 0 eccentricity. Do not pass a ''center'' param.');
                    end
                    vararginForConeMosaic{numel(vararginForConeMosaic)+1} = varargin{k};
                    vararginForConeMosaic{numel(vararginForConeMosaic)+1} = varargin{k+1};
                end
            end
            obj = obj@coneMosaic(vararginForConeMosaic{:});
    
            % parse input
            p = inputParser;
            p.addRequired('resamplingFactor', @isnumeric);
            p.addParameter('fovDegs', 0.25, @(x)(isnumeric(x)&&((numel(x)==1)||(numel(x)==2))));
            p.addParameter('eccBasedConeDensity', false, @islogical);
            p.addParameter('sConeMinDistanceFactor', 3.0, @isnumeric);
            p.addParameter('sConeFreeRadiusMicrons', 45, @isnumeric);
            p.addParameter('customInnerSegmentDiameter', [], @isnumeric);
            p.addParameter('customLambda', [], @isnumeric);
            p.addParameter('rotationDegs', 0, @isnumeric);
            p.addParameter('latticeAdjustmentPositionalToleranceF', 0.01, @isnumeric);
            p.addParameter('latticeAdjustmentDelaunayToleranceF', 0.001, @isnumeric);
            p.addParameter('saveLatticeAdjustmentProgression', false, @islogical);
            p.addParameter('marginF', 1.75, @(x)(isnumeric(x)&&(x>0.0)));
            p.parse(upSampleFactor, vararginForConeHexMosaic{:});
            
            % Set input params
            obj.resamplingFactor = p.Results.resamplingFactor;
            obj.eccBasedConeDensity = p.Results.eccBasedConeDensity;
            obj.sConeMinDistanceFactor = p.Results.sConeMinDistanceFactor;
            obj.sConeFreeRadiusMicrons = p.Results.sConeFreeRadiusMicrons;
            obj.customLambda = p.Results.customLambda;
            obj.customInnerSegmentDiameter = p.Results.customInnerSegmentDiameter;
            obj.rotationDegs = p.Results.rotationDegs;
            obj.saveLatticeAdjustmentProgression = p.Results.saveLatticeAdjustmentProgression;
            obj.latticeAdjustmentDelaunayToleranceF = p.Results.latticeAdjustmentDelaunayToleranceF;
            obj.latticeAdjustmentPositionalToleranceF = p.Results.latticeAdjustmentPositionalToleranceF;
            
            % Set FOV of the underlying rect mosaic
            if (numel(p.Results.fovDegs) == 1)
                obj.setSizeToFOV(p.Results.fovDegs(1)*[1 1]);
            else
                obj.setSizeToFOV(p.Results.fovDegs);
            end
            
            if (isempty(p.Results.marginF))
                if (max(p.Results.fovDegs) >= 1.0)
                    obj.marginF = 1.5;
                else
                    obj.marginF = 1.8;
                end
            else
                obj.marginF = p.Results.marginF;
            end
            
            % Get a copy of the original coneLocs
            obj.saveOriginalResState();
            
            % Set custom pigment light collecting dimensions
            if (~isempty(obj.customInnerSegmentDiameter)) 
                maxInnerSegmentDiameter = 1e6 * diameterForCircularApertureFromWidthForSquareAperture(obj.pigment.pdWidth);
                if (obj.eccBasedConeDensity) && (obj.customInnerSegmentDiameter>maxInnerSegmentDiameter)
                   error('The custom inner segment diameter (%2.4f) is > max inner segment diameter (%2.4f) necessary to keep the default cone density. Either set ''eccBasedConeDensity'' to false, or decrease ''customInnerSegmentDiameter''.', obj.customInnerSegmentDiameter, maxInnerSegmentDiameter);
                end
                obj.pigment.pdWidth = 1e-6 * sizeForSquareApertureFromDiameterForCircularAperture(obj.customInnerSegmentDiameter);
                obj.pigment.pdHeight = obj.pigment.pdWidth;
            end
            
            % Set the pigment geometric dimensions
            if (~isempty(obj.customLambda)) 
                maxSpacing = 1e6 * diameterForCircularApertureFromWidthForSquareAperture(obj.pigment.width);
                if (obj.eccBasedConeDensity) && (obj.customLambda > maxSpacing)
                    error('The custom lambda (%2.4f) is > max separation (%2.4f) in order to keep the default cone density. Either set ''eccBasedConeDensity'' to false, or decrease ''customLambda''.', obj.customLambda, maxSpacing);
                end
                obj.pigment.width = 1e-6 * obj.customLambda;
                obj.pigment.height = obj.pigment.width;
            end
            
            % Generate sampled hex grid
            obj.resampleGrid(obj.resamplingFactor);
            
            if ((~isempty(obj.sConeMinDistanceFactor)) || (~isempty(obj.sConeFreeRadiusMicrons)))
                % Make s-cone lattice semi-regular, and/or add an s-cone free region.
                obj.reassignConeIdentities(...
                    'sConeMinDistanceFactor', obj.sConeMinDistanceFactor, ...   
                    'sConeFreeRadiusMicrons', obj.sConeFreeRadiusMicrons);    
            end
        end
        
        % Visualize different aspects of the hex grid
        hFig = visualizeGrid(obj, varargin);
        
        % Method to compute the cone density of @coneMosaicHex
        [densityMap, densityMapSupportX, densityMapSupportY] = computeDensityMap(obj, computeConeDensityMap)
        
        % Reshape a full 3D hex activation map (coneRows x coneCols x time] to a 2D map (non-null cones x time)
        hex2Dmap = reshapeHex3DmapToHex2Dmap(obj, hex3Dmap);
        
        % Reshape a 2D map (non-null cones x time) to the full 3D hex activation map (coneRows x coneCols x time)
        hex3Dmap = reshapeHex2DmapToHex3Dmap(obj, hex2Dmap);
        
        % Compute activation images for the hex mosaic (all cones +  LMS submosaics)
        [activationImage, activationImageLMScone, imageXaxis, imageYaxis] = computeActivationDensityMap(obj, activation);
        
        % Visualize activation map images for the hex mosaic (all cones +  LMS submosaics)
        hFig = visualizeActivationMaps(obj, activation, varargin);
        
        % Render (in the passed axesHandle) an activation map for the hex mosaic
        renderActivationMap(obj, axesHandle, activation, varargin);
        
        % Visualize iterative adjustment of the cone lattice 
        hFig = plotMosaicProgression(obj, varargin);
        
        % Print various infos about the cone mosaic
        displayInfo(obj);
    end % Public methods
    
    methods (Access = private)
        % Private methods
        saveOriginalResState(obj);
        restoreOriginalResState(obj);
        
        % Change the FOV of the mosaic
        setSizeToFOVForHexMosaic(obj,fov);
        
        % Change the cone identities according to arguments passed in varargin
        reassignConeIdentities(obj, varargin);
        
        % Sample the original rectangular mosaic using a hex grid sampled at the passed resamplingFactor
        resampleGrid(obj, resamplingFactor);
    end % Private methods
    
    methods (Static)
        renderPatchArray(axesHandle, pixelOutline, xCoords, yCoords, edgeColor, faceColor, lineStyle);
        renderHexMesh(axesHandle, xHex, yHex, meshEdgeColor, meshFaceColor, meshFaceAlpha, meshEdgeAlpha, lineStyle);
    end % Static methods
    
end

