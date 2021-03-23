classdef cMosaic < handle
    % Create an optimized cone mosaic object
    %
    
    % Syntax:
    %   theConeMosaic = cMosaic;
    %   theConeMosaic = cMosaic('pigment', pp);
    %
    % Description:
    %    An object of the @cMosaic class defines parameters of the
    %    retinal cone mosaic, and enables computation of cone
    %    isomerizations and photocurrent response in an array of cones.
    %
    %    theConeMosaic =  cMosaic(..., 'PARAM1', val1, 'PARAM2', val2, ...)
    %    creates the cone mosaic object. Optional parameter name/value
    %    pairs are listed below.
    %
    %    ...
    %
    % Inputs:
    %    None required.
    %
    % Outputs:
    %    cMosaic           The created coneMosaic object.
    %
    % Optional key/value pairs:
    %    'name'             - String. Mosaic name. Default 'cone mosaic'.
    %
    % See Also:
    %    cPhotoPigment
    %
    % History:
    %    December 2020  NPC  Wrote it
    
    properties (Constant)
        % Cone types
        LCONE_ID = 1;
        MCONE_ID = 2;
        SCONE_ID = 3;
        KCONE_ID = 4;
        
        % Cone aperture (light collecting disk) has a diameter that 
        % is 0.8 x cone spacing. We chose this to be consistent with
        % coneMosaicHex, .e.g.: 
        % c = coneMosaicHex(3);
        % c.pigment.pdWidth/c.pigment.width, which gives 0.7899
        coneApertureToDiameterRatio = 0.79;
    end
    
    % Public properties
    properties  (GetAccess=public, SetAccess=public)
        % Name of the mosaic
        name;
        
        % Cone photopigment for the mosaic
        pigment;
        
        % Macular pigment object for mosaic
        macular;
        
        % Boolean indicating whether to account for the increased photon
        % capture due to the increase in cone aperture with eccentricity
        eccVaryingConeAperture;
        
        % Boolean indicating whether to blur by eccentricity-varying cone apertures 
        % (true) of the median cone aperture (false). This quantization of
        % the eccentricity-varying cone apertures will depend on the pixel
        % size of the input optical image. The smallest the pixel size, the
        % finer the quantization.
        eccVaryingConeBlur;

        % Boolean indicating whether to account for the decrease in outer segment 
        % cone length with increasing eccentricity.
        eccVaryingOuterSegmentLength;
        
        % Boolean indicating whether to account for the decrease in macular pigment 
        % density with increasing eccentricity.
        eccVaryingMacularPigmentDensity;
        
        % Boolean indicating whether to account for the decrease in macular pigment 
        % density with increasing eccentricity in the presence of eye
        % movements
        eccVaryingMacularPigmentDensityDynamic;
        
        % 3- or 4-element vector containing relative cone densities
        coneDensities;
        
        % Radius of S-cone free area around (0,0)
        tritanopicRadiusDegs;
        
        % Poisson noise flag for cone excitations 
        noiseFlag;
        
        % Integration time (seconds)
        integrationTime;
        
        % Random seed
        randomSeed;
        
        % Color for cone rendering in mosaic visualization
        lConeColor = [1 0.2 0.3];
        mConeColor = [0.2 1 0.4];
        sConeColor = [0.3 0.1 1];
        kConeColor = [0.9 0.9 0.2];
    end
    
    % Read-only properties
    properties (GetAccess=public, SetAccess=private)
        % [n x 2] matrix of RGC positions, in microns
        coneRFpositionsMicrons;
        
        % [n x 2] matrix of RGC positions, in degrees
        coneRFpositionsDegs;
        
        % [n x 1] vector of cone spacings, in degrees
        coneRFspacingsDegs;
        
        % [n x 1] vector of cone spacings, in microns
        coneRFspacingsMicrons;
        
        % [n x 1] vector of cone types
        coneTypes;
        
        % indices of L-cones
        lConeIndices;
        % indices of M-cones
        mConeIndices;
        % indices of S-cones
        sConeIndices;
        % indices of K-cones
        kConeIndices;
        
        % Eccentricity [x,y] of the center of the cone mosaic 
        eccentricityDegs;
        
        % Size of the cone mosaic [width, height]
        sizeDegs;
        
        % Either 'right eye' or 'left eye'
        whichEye;
        
        % Photon absorption attenuation factors to account for decrease in 
        % outer segment length with eccentricity
        outerSegmentLengthEccVariationAttenuationFactors = [];
        
        % Fixational eye movement object for the mosaic
        fixEMobj = [];
        
        % Blur diameter for the different zones of cones
        blurApertureDiameterMicronsZones = [];
        
        % User-settable microns per degree
        micronsPerDegreeApproximation = [];
    end
    
    % Dependent properties
    properties (Dependent)
        % wave - Wavelength support
        %   wave depends on wavelength support of the pigment property
        wave;
        
        % qe - Cone quantal efficiency (absorptance)
        %   The quantal efficiency is the product of the cone
        %   photopigment absorptance times the macular pigment
        %   transmittance times the cone photopigment peak efficiency.
        %
        %   The quantal efficiency referred all the way to the cornea and
        %   used in the plot function, also includes the lens
        %   transmittance. The lens transmittance is specified by the LENS
        %   object in the oi representation.
        %
        %   The QE depends on the PIGMENT and MACULAR properties.
        qe;
        
        % average microns per degree
        micronsPerDegree;
        
        % mosaic size in microns
        sizeMicrons;
    end
    
    
    % Private properties
    properties (GetAccess=private, SetAccess=private)
        % Size (in degs) of source lattice from which to crop positions for
        % the desired eccentricity
        sourceLatticeSizeDegs = 45;
        
        % Struct containing information related to the current cone partition
        % into zones based on their aperture size
        cachedConePartition = [];
        
        % Struct containing information related to the current
        % macularPigmentDensityBoostFactors
        cachedMacularPigmentDensityBoostFactors = [];
        
        % OSLength attenuation factors (only used when importing coneData)
        importedOSLengthAttenuationFactors = [];
        
        % cone aperture blur (only used when importing coneData)
        importedBlurDiameterMicrons = [];
        
        % Whether to use parfor in computations
        useParfor = true;
        
        % Min and max cone positions
        minRFpositionMicrons;
        maxRFpositionMicrons;
    end
    
    % Public methods
    methods
        
        % Constructor
        function obj = cMosaic(varargin)
            % Parse input
            p = inputParser;
            p.addParameter('name', 'cone mosaic', @ischar);
            p.addParameter('wave', 400:10:700, @isnumeric);
            p.addParameter('pigment', cPhotoPigment(), @(x) isa(x, 'cPhotoPigment'));
            p.addParameter('macular', Macular(), @(x)isa(x, 'Macular'));
            p.addParameter('coneData', [], @(x)(isempty(x) || (isstruct(x))));
            p.addParameter('eccentricityDegs', [0 0], @(x)(isnumeric(x) && (numel(x) == 2)));
            p.addParameter('sizeDegs', [0.4 0.4], @(x)(isnumeric(x) && (numel(x) == 2)));
            p.addParameter('computeMeshFromScratch', false, @islogical);
            p.addParameter('visualizeMeshConvergence', false, @islogical);
            p.addParameter('exportMeshConvergenceHistoryToFile', false, @islogical);
            p.addParameter('maxMeshIterations', 100, @(x)(isempty(x) || isscalar(x)));
            p.addParameter('whichEye', 'right eye', @(x)(ischar(x) && (ismember(x, {'left eye', 'right eye'}))));
            p.addParameter('micronsPerDegree', [], @(x)(isempty(x) || (isscalar(x))));
            p.addParameter('eccVaryingConeAperture', true, @islogical);
            p.addParameter('eccVaryingConeBlur', false, @islogical);
            p.addParameter('eccVaryingOuterSegmentLength', true, @islogical);
            p.addParameter('eccVaryingMacularPigmentDensity', true, @islogical);
            p.addParameter('eccVaryingMacularPigmentDensityDynamic', false, @islogical);
            p.addParameter('coneDensities', [0.6 0.3 0.1 0.0], @(x)(isnumeric(x) && ((numel(x) == 3)||(numel(x)==4))));
            p.addParameter('tritanopicRadiusDegs', 0.15, @isscalar);
            p.addParameter('noiseFlag', 'random', @(x)(ischar(x) && (ismember(x, {'random', 'frozen', 'none'}))));
            p.addParameter('randomSeed', [], @isscalar);
            p.addParameter('integrationTime', 5/1000, @isscalar);
            p.addParameter('useParfor', true, @islogical);
            
            p.parse(varargin{:});
            
            obj.name = p.Results.name;
            obj.macular = p.Results.macular;
            obj.pigment = p.Results.pigment;
            obj.wave = p.Results.wave;
            
            obj.eccentricityDegs = p.Results.eccentricityDegs;
            obj.sizeDegs = p.Results.sizeDegs;
            obj.whichEye = p.Results.whichEye;

            obj.micronsPerDegreeApproximation = p.Results.micronsPerDegree;
            obj.eccVaryingConeAperture = p.Results.eccVaryingConeAperture;
            obj.eccVaryingOuterSegmentLength = p.Results.eccVaryingOuterSegmentLength;
            obj.eccVaryingMacularPigmentDensity = p.Results.eccVaryingMacularPigmentDensity;
            
            % New capabilities for more realistic but slower simulation
            obj.eccVaryingConeBlur = p.Results.eccVaryingConeBlur;
            obj.eccVaryingMacularPigmentDensityDynamic  = p.Results.eccVaryingMacularPigmentDensityDynamic;
            
            obj.coneDensities = p.Results.coneDensities;
            obj.tritanopicRadiusDegs = p.Results.tritanopicRadiusDegs;
            
            obj.noiseFlag = p.Results.noiseFlag;
            obj.randomSeed = p.Results.randomSeed;
            obj.integrationTime = p.Results.integrationTime;
            
            % Parallel computations
            obj.useParfor = p.Results.useParfor;
            
            % Assert that we have appropriate pigment if we have more than
            % 3 cone types
            if (numel(obj.coneDensities)>3) && (any(obj.coneDensities(4:end)>0.0))
                assert(numel(obj.coneDensities) == size(obj.pigment.absorptance,3), ...
                    sprintf('cPhotoPigment is not initialized for %d types of cones', numel(obj.coneDensities)));
            end
            
            % These listeners make sure the wavelength support
            % in obj.pigment and obj.macular match the wave property
            addlistener(obj.pigment, 'wave', 'PostSet', @obj.matchWaveInAttachedMacular);
            addlistener(obj.macular, 'wave', 'PostSet', @obj.matchWaveInAttachedPhotopigment);
            
            if (isempty(p.Results.coneData))
                if (p.Results.computeMeshFromScratch)
                    % Re-generate lattice
                    visualizeMeshConvergence = p.Results.visualizeMeshConvergence; 
                    exportMeshConvergenceHistoryToFile = p.Results.exportMeshConvergenceHistoryToFile;
                    obj.regenerateConePositions(p.Results.maxMeshIterations,  ...
                        visualizeMeshConvergence, exportMeshConvergenceHistoryToFile);
                else
                    % Import positions by cropping a large pre-computed patch
                    obj.initializeConePositions();
                end
                
                % Set random seed
                if (isempty(obj.randomSeed))
                    rng('shuffle');
                else
                    rng(obj.randomSeed);
                end

                % Assign cone types
                obj.assignConeTypes();
            else
                % Load cone positions and cone types from passed struct
                obj.importExternalConeData(p.Results.coneData);

                % Set random seed
                if (isempty(obj.randomSeed))
                    rng('shuffle');
                else
                    rng(obj.randomSeed);
                end
            end
            
            % Save min and max cone position
            obj.minRFpositionMicrons = squeeze(min(obj.coneRFpositionsMicrons,[],1));
            obj.maxRFpositionMicrons = squeeze(max(obj.coneRFpositionsMicrons,[],1));
            
            % Compute photon absorption attenuation factors to account for
            % the decrease in outer segment legth with ecc.
            obj.computeOuterSegmentLengthEccVariationAttenuationFactors('useParfor', obj.useParfor);
        end
        
        % Method to assign cone types
        assignConeTypes(obj, varargin);
        
        % Method to visualize the cone mosaic
        visualize(obj, varargin);
        
        % Method to generate an eye movement sequence
        emGenSequence(obj, durationSeconds, varargin);
        
        % Method to generate an ensemble of OIs for the mosaic
        [oiEnsemble, psfEnsemble] = oiEnsembleGenerate(obj, oiSamplingGridDegs, varargin);
        
        % Method to compute the mosaic response
        [absorptionsCount, noisyAbsorptionInstances, ...
            photoCurrents, photoCurrentInstances, responseTemporalSupport]  = compute(obj, oi, varargin);
        
        % Method for suggesting a minimal pixel size (in degrees) for the input scene
        % depending on the aperture size of the mosaic's cones and whether
        % the computation is done using ecc-dependent blur mode
        scenePixelSizeDegs = suggestedScenePixelSizeDegs(obj, eccVaryingConeBlur);
        
        % Getter/Setter methods for dependent variables
        % QE
        function val = get.qe(obj)
            % Compute the quantal efficiency of the cones
            val = bsxfun(@times, obj.pigment.absorptance, obj.macular.transmittance) * ...
                  diag(obj.pigment.peakEfficiency);
        end
        
        % WAVE - getter
        function val = get.wave(obj)
            val = obj.pigment.wave;
        end
        % WAVE - setter
        function set.wave(obj, val)
            obj.pigment.wave = val(:);
            obj.macular.wave = val(:);
        end
        
        % MICRONSPERDEGREE
        function val = get.micronsPerDegree(obj)
            val = (max(obj.coneRFpositionsMicrons,[],1) - min(obj.coneRFpositionsMicrons,[],1)) / ...
                  (max(obj.coneRFpositionsDegs,[],1)    - min(obj.coneRFpositionsDegs,[],1));
        end
        
        % SIZEMICRONS
        function val = get.sizeMicrons(obj)
            xyMax = max(obj.coneRFpositionsMicrons,[],1);
            xyMin = min(obj.coneRFpositionsMicrons,[],1);
            val = xyMax - xyMin;
        end
    end
    
    
    methods (Access=private)
        
        % Initialize cone positions by importing them from a large previously-computed mesh
        initializeConePositions(obj);
        
        % Initialize cone positions by regenerating a new mesh. Can be slow.
        regenerateConePositions(obj, maxIterations, visualizeConvergence, exportHistoryToFile);
        
        % Method to partition cones into zones based on cone
        % aperture size and current optical image resolution
        [coneApertureDiameterMicronsZoneBands, ...     % the median cone aperture in this zone band
         coneIndicesInZones] =  ...                    % the IDs of cones in this zone band
            coneZonesFromApertureSizeAndOIresolution(obj, coneApertureDiametersMicrons, oiResMicrons);
        
        % Method to compute boost factors for correction 
        macularPigmentDensityBoostFactors = computeMPBoostFactors(obj, oiPositionsDegs, emPositionDegs, oiWave, oiSize, oiResMicrons)
        
        % Compute absorption rate
        absorptionsRates = computeAbsorptionRate(obj, currentEMposMicrons, ...
            oiPositionsMicrons, absorptionsDensityImage,  ...
            oiResMicrons, coneApertureDiametersMicrons, ...
            coneIndicesInZones);

        % Called when the wave property in the attached photopigment object is changed
        function matchWaveInAttachedMacular(obj, src,~)
            % Set the wave property in the attached macular
            obj.macular.wave = obj.pigment.wave;
        end
        
        % Called when the wave property in the attached macular object is changed
        function matchWaveInAttachedPhotopigment(obj, src,~)
            % Set the wave property in the attached cPhotoPigment
            obj.pigment.wave = obj.macular.wave;
        end
    end
    
    methods (Static)
        % Static method to generate cone aprture blur kernel
        apertureKernel = generateApertureKernel(coneApertureDiameterMicrons, oiResMicrons);
        % Static method to generate noisy absorption response instances
        % from the mean absorption responses
        noisyAbsorptionInstances = noisyInstances(meanAbsorptions, varargin);
    end
end

