classdef mRGCMosaic < handle
% Create a midget RGC mosaic connected to a @cMosaic

    % Public properties
    properties  (GetAccess=public, SetAccess=public)
        % Name of the mosaic
        name;
        
        % Poisson noise flag for cone excitations 
        noiseFlag;
        
        % Random seed
        randomSeed;
    end
    
    % Read-only properties
    properties (GetAccess=public, SetAccess=private)
        % [n x 2] matrix of RGC positions, in microns
        rgcRFpositionsMicrons;
        
        % [n x 2] matrix of RGC positions, in degrees
        rgcRFpositionsDegs;
        
        % [n x 1] vector of RGC spacings, in degrees
        rgcRFspacingsDegs;
        
        % [n x 1] vector of rRGC spacings, in microns
        rgcRFspacingsMicrons;
        
        % Eccentricity [x,y] of the center of the RGC mosaic (degs)
        eccentricityDegs;
        
        % Eccentricity [x,y] of the center of the RGC mosaic (microns)
        eccentricityMicrons; 
        
        % Size of the RGC mosaic [width, height]
        sizeDegs;
        
        % Either 'right eye' or 'left eye'
        whichEye;
        
        % Generated input cone mosaic
        inputConeMosaic;
        
        % Passed from input cone mosaic
        micronsPerDegreeApproximation;
    end
    
    % Private properties
    properties (GetAccess=private, SetAccess=private)
        % Size (in degs) of source lattice from which to crop positions for
        % the desired eccentricity
        sourceLatticeSizeDegs = 58;
        
        % Whether to use parfor in computations
        useParfor = true;
        
        % Min and max cone positions
        minRFpositionMicrons;
        maxRFpositionMicrons;
        
        % Sparse [conesNum x rgcsNum] sparse  connectivity matrix 
        % of cone -> rgc. 
        % To find which cones are connected to a target RGC:
        %  connectivityVector = full(squeeze(obj.coneConnectivityMatrix(:, targetRGC)));
        %  inputConeIDs = find(connectivityVector > 0.01);
        coneConnectivityMatrix;
    end
    
    % Public methods
    methods
        % Constructor
        function obj = mRGCMosaic(varargin)
            % Parse input
            p = inputParser;
            p.addParameter('name', 'mRGC mosaic', @ischar);
            p.addParameter('whichEye', 'left eye', @(x)(ischar(x) && (ismember(x, {'left eye', 'right eye'}))));
            p.addParameter('eccentricityDegs', [0 0], @(x)(isnumeric(x) && (numel(x) == 2)));
            p.addParameter('sizeDegs', [0.4 0.4], @(x)(isnumeric(x) && (numel(x) == 2)));
            p.addParameter('computeMeshFromScratch', false, @islogical);
            p.addParameter('maxMeshIterations', 100, @(x)(isempty(x) || isscalar(x)));
            p.addParameter('visualizeMeshConvergence', false, @islogical);
            p.addParameter('exportMeshConvergenceHistoryToFile', false, @islogical);
            p.addParameter('inputConeMosaic', [], @(x)(isempty(x) || isa(x, 'cMosaic')));
            p.addParameter('noiseFlag', 'random', @(x)(ischar(x) && (ismember(x, {'random', 'frozen', 'none'}))));
            p.addParameter('randomSeed', [], @isscalar);
            p.parse(varargin{:});
            
            obj.name = p.Results.name;
            obj.eccentricityDegs = p.Results.eccentricityDegs;
            obj.sizeDegs = p.Results.sizeDegs;
            obj.whichEye = p.Results.whichEye;
            
            obj.noiseFlag = p.Results.noiseFlag;
            obj.randomSeed = p.Results.randomSeed;
            
            % Set the inputConeMosaic
            if (isempty(p.Results.inputConeMosaic))
                % Generate input cone mosaic (make a bit larger so we have
                % enough cones for the surrounds of the most eccentric RGCs)
                maxEccentricityDegs = sqrt(sum((obj.eccentricityDegs + 0.5*obj.sizeDegs).^2,2));
                extraDegsForRGCSurround = 2.0 * RGCmodels.CronerKaplan.constants.surroundCharacteristicRadiusFromFitToPandMcells(maxEccentricityDegs);

                obj.inputConeMosaic = cMosaic(...
                    'whichEye', obj.whichEye, ...
                    'sizeDegs', obj.sizeDegs + 2 * extraDegsForRGCSurround, ...
                    'eccentricityDegs', obj.eccentricityDegs, ...
                    'opticalImagePositionDegs', 'mosaic-centered');
            else
                obj.inputConeMosaic = p.Results.inputConeMosaic;
            end
            
            % Set the microns per deg approximation
            obj.micronsPerDegreeApproximation = obj.inputConeMosaic.micronsPerDegreeApproximation;
            
            if (p.Results.computeMeshFromScratch)
                % Re-generate lattice
                obj.regenerateRFPositions(...
                        p.Results.maxMeshIterations,  ...
                        p.Results.visualizeMeshConvergence, ...
                        p.Results.exportMeshConvergenceHistoryToFile);
            else
                % Import positions by cropping a large pre-computed patch
                obj.initializeRFPositions();
            end
            
            % Remove RFs within the optic disk
            obj.removeRFsWithinOpticNerveHead();
                
            % Set random seed
            if (isempty(obj.randomSeed))
               rng('shuffle');
            else
               rng(obj.randomSeed);
            end

            % Compute min and max cone position
            obj.minRFpositionMicrons = squeeze(min(obj.rgcRFpositionsMicrons,[],1));
            obj.maxRFpositionMicrons = squeeze(max(obj.rgcRFpositionsMicrons,[],1));
            
            % Compute ecc in microns
            if (isempty(obj.minRFpositionMicrons))
                % No rgcs, set value differently
                if (~isempty(obj.micronsPerDegreeApproximation))
                    obj.eccentricityMicrons = obj.eccentricityDegs * obj.micronsPerDegreeApproximation;
                else
                    obj.eccentricityMicrons = obj.eccentricityDegs *  300;
                end
            else
                obj.eccentricityMicrons = 0.5*(obj.minRFpositionMicrons + obj.maxRFpositionMicrons);
            end
            
            % Wire mRGCs RF centers to cones of the input cone mosaic
            visualizeAlignment = false; 
            visualizeWiringStages = false;
            obj.wireRFcentersToInputCones(visualizeAlignment, visualizeWiringStages);
        end
        
        % Method to return indices of RFs within an ROI
        rfIndices = indicesOfRFsWithinROI(obj, roi);
        
        % Method to return the cone-to-RGC density ratio for all RGCs
        densityRatios = coneToRGCDensityRatios(obj);
        
        % Method to visualize the mRGCmosaic
        visualize(obj, varargin);
        
        % Method to visualize how the cone mosaic is tesselated by the RGC RF centers
        visualizeConeMosaicTesselation(obj, ...
            coneRFpositions, coneRFspacings, ...
            rgcRFpositions, rgcRFspacings, ...
            domain, varargin);
        
        % Method to report the distribution of cones / RF center
        connectivityStats(obj, figNo);
    end
    
    methods (Access=private)
        % Initialize RF positions by importing them from a large previously-computed mesh
        initializeRFPositions(obj);
        
        % Initialize RF positions by regenerating a new mesh. Can be slow.
        regenerateRFPositions(obj, maxIterations, visualizeConvergence, exportHistoryToFile);
        
        % Remove RFs located within the optic disk
        removeRFsWithinOpticNerveHead(obj);
        
        % Method to wire RF centers to cones of the input cone mosaic
        wireRFcentersToInputCones(obj, visualizeAlignment, visualizeWiringStages);
        
        % Method to align RGC RFs to cones in the central retina.
        % Called by obj.wireRFcenterToInputCones()
        alignRGCsWithCones(obj, coneRFpositionsMicrons, coneRFpositionsDegs, visualizeAlignment)
    
        % Method to connect cones to RGC RF centers
        % Called by obj.wireRFcenterToInputCones()
        connectConesToRGCcenters(obj, idxConesInsideRGCmosaic, visualizeConnection)
    
    end
    
end

