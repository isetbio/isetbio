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
    % Inputs:
    %    None required.
    %
    % Outputs:
    %    cMosaic           The created cMosaic object.
    %
    % Optional key/value pairs:
    %    'name'                             - String. Mosaic name. Default 'cone mosaic'.
    %    'wave'                             - Vector. Wavalength support. Default is: 400:10:700.
    %    'pigment'                          - @cPhotoPigment. Custom cPhotoPigment object.
    %    'macular'                          - @Macular(). Custom Macular object.
    %    'coneData',                        - Struct. Struct with custom cone data, usually exported from a @coneMosaicHex.
    %    'eccentricityDegs'                 - [x,y]. Center of the mosaic, in degrees. Default: [0 0]
    %    'sizeDegs'                         - [sx, sy]. Width and height of the mosaic, in degrees. Default: [0.4 0.4]
    %    'whichEye'                         - Which eye. Valid options are: {'left eye', 'right eye'}. Default: 'right eye'.
    %    'computeMeshFromScratch'           - Logical, indicating whether to recompute the mesh from scratch. Default: false
    %    'visualizeMeshConvergence'         - Logical, indicating whether to visualize the convergence of the mesh. Default: false 
    %    'exportMeshConvergenceHistoryToFile' - Logical, indicating whether to save the convergence of the mesh. Default:false
    %    'maxMeshIterations'                - Scalar. Max number of iterations for the mesh generation. Default: 100.
    %    'customRFspacingFunction'          - Empty [], of a handle to a function that returns the RF spacing in microns for an array of RF positions in microns - used only when generating a mesh from scratch 
    %    'customDegsToMMsConversionFunction'- Empty [], of a handle to a function that returns eccentricity in degs of visual angle for eccentricities specified in to retinal mms - used only when generating a mesh from scratch
    %    'customMMsToDegsConversionFunction'- Empty [], of a handle to a function that returns eccentricity in retinal mms for eccentricities specified in to degrees of visual angle   - used only when generating a mesh from scratch        
    %    'micronsPerDegree'                 - Scalar. A custom retinal magnification factor, microns/deg
    %    'eccVaryingConeAperture'           - Logical, indicating whether to allow the cone aperture (light collecting area) to vary with eccentricity. Default: true 
    %    'eccVaryingConeBlur'               - Logical, indicating whether to allow the cone aperture (spatial blur) to vary with eccentricity.  Default: false 
    %    'eccVaryingOuterSegmentLength'     - Logical, indicating whether to allow the outer segment (light collecting area) to vary with eccentricity. Default: true
    %    'eccVaryingMacularPigmentDensity'  - Logical, indicating whether to allow the macular pigment density to vary with eccentricity. Default: true
    %    'eccVaryingMacularPigmentDensityDynamic',- Logical, indicating whether to correct the macular pigment density in the presence of eye movements. Default: false
    %    'coneDensities'                    - Vector, with 3 or 4 densities for L-, M-, S-, and possibly a 4th cone. Default: [0.6 0.3 0.1 0.0]
    %    'tritanopicRadiusDegs'             - Scalar. Radius of S-cone free region. Default: 0.15 degs
    %    'noiseFlag'                        - String. Indicating whether to compute with random noise, frozen noise or no-noise. Valid options are {'random', 'frozen', 'none'}
    %    'randomSeed'                       - Random seed for the noise.
    %    'integrationTime'                  - Scalar. Integration time in seconds. Default is: 5/1000
    %    'useParfor'                        - Logical, indicating whether to use parfor. Default: true.
    %
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
        % is 0.79 x cone spacing. We chose this to be consistent with
        % coneMosaicHex, .e.g.: 
        % c = coneMosaicHex(3);
        % c.pigment.pdWidth/c.pigment.width, which gives 0.7899
        coneApertureToDiameterRatio = 0.79;
        
        % OD size is from "The Size and Shape of the Optic Disc in Normal Human Eyes"
        %                  Harry A. Quigley, MD; Andrew E. Brown; John D. Morrison, MD; et al Stephen M. Drance, MD
        %                  Arch Ophthalmol. 1990;108(1):51-57. doi:10.1001/archopht.1990.01070030057028
        %
        % OD location is from "Determination of the Location of the Fovea on the Fundus".
        %                  Invest. Ophthalmol. Vis. Sci. 2004;45(9):3257-3258. doi: https://doi.org/10.1167/iovs.03-1157.
        % These are both specified in retinal coordinates
        opticDisk = struct(...
            'centerDegs', [15.5, 1.5], ...
        	'horizontalDiameterMM', 1.77, ...
        	'verticalDiameterMM', 1.88, ...
            'rotationDegs', 6.0);
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
        
        % Poisson noise flag for cone excitations 
        noiseFlag;
        
        % Integration time (seconds)
        integrationTime;
        
        % Random seed
        randomSeed;
        
        % Default optical image position degs. This can be overriden by
        % passing a ('opticalImagePositionDegs', something) key/value pair during the cMosaic.compute() call.
        opticalImagePositionDegs;
        
        % Color for cone rendering in mosaic visualization
        lConeColor = [1 0.1 0.5];
        mConeColor = [0.1 1 0.5];
        sConeColor = [0.6 0.1 1];
        kConeColor = [0.9 0.9 0.2];
    end
    
    
    properties (SetObservable, AbortSet)
        % These properties trigger the assignConeTypes() method
        % 3- or 4-element vector containing relative cone densities
        coneDensities;
        
        % Radius of S-cone free area around (0,0)
        tritanopicRadiusDegs;
    end
    
    % Read-only properties
    properties (GetAccess=public, SetAccess=private)
        % [n x 2] matrix of cone positions, in microns
        coneRFpositionsMicrons;
        
        % [n x 2] matrix of cone positions, in degrees
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
        
        % Eccentricity [x,y] of the center of the cone mosaic (degs)
        eccentricityDegs;
        
        % Eccentricity [x,y] of the center of the cone mosaic (microns)
        eccentricityMicrons; 
        
        % Size of the cone mosaic [width, height]
        sizeDegs;
        
        % Either 'right eye' or 'left eye'
        whichEye;
        
        % Fixational eye movement object for the mosaic
        fixEMobj = [];
        
        % Blur diameter for the different zones of cones
        blurApertureDiameterMicronsZones = [];
        
        % User-settable microns per degree
        micronsPerDegreeApproximation = [];
        
        % Achieved cone densities
        achievedConeDensities;
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
        sourceLatticeSizeDegs = 58;
        
        % Struct containing information related to the current cone partition
        % into zones based on their aperture size
        cachedConePartition = [];
        
        % Struct containing information related to the current
        % macularPigmentDensityBoostFactors
        cachedMacularPigmentDensityBoostFactors = [];
        
        % OSLength attenuation factors (only used when importing coneData)
        importedOSLengthAttenuationFactors = [];
        
        % Photon absorption attenuation factors to account for decrease in 
        % outer segment length with eccentricity
        outerSegmentLengthEccVariationAttenuationFactors = [];
        
        % cone aperture blur (only used when importing coneData)
        importedBlurDiameterMicrons = [];
        
        % Whether to use parfor in computations
        useParfor = true;
        
        % Min and max cone positions
        minRFpositionMicrons;
        maxRFpositionMicrons;
        
        % Full absorptions density map
        absorptionsDensityFullMap;
        absorptionsDensitySpatialSupportMicrons;
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
            p.addParameter('whichEye', 'right eye', @(x)(ischar(x) && (ismember(x, {'left eye', 'right eye'}))));
            p.addParameter('computeMeshFromScratch', false, @islogical);
            p.addParameter('customMinRFspacing', [], @(x) (isempty(x) || isscalar(x)));
            p.addParameter('customRFspacingFunction', [], @(x) (isempty(x) || isa(x,'function_handle')));
            p.addParameter('customDegsToMMsConversionFunction', [], @(x) (isempty(x) || isa(x,'function_handle')));
            p.addParameter('customMMsToDegsConversionFunction', [], @(x) (isempty(x) || isa(x,'function_handle')));
            p.addParameter('visualizeMeshConvergence', false, @islogical);
            p.addParameter('exportMeshConvergenceHistoryToFile', false, @islogical);
            p.addParameter('maxMeshIterations', 100, @(x)(isempty(x) || isscalar(x)));
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
            p.addParameter('opticalImagePositionDegs', 'mosaic-centered', @(x)(ischar(x) || (isnumeric(x)&&numel(x)==2)));
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
            obj.opticalImagePositionDegs = p.Results.opticalImagePositionDegs;
            
            % Parallel computations
            obj.useParfor = p.Results.useParfor;
            
            % Custom mesh generation function 
            customMinRFspacing = p.Results.customMinRFspacing;
            customRFspacingFunction = p.Results.customRFspacingFunction;
            customDegsToMMsConversionFunction = p.Results.customDegsToMMsConversionFunction;
            customMMsToDegsConversionFunction = p.Results.customMMsToDegsConversionFunction;
            
            % Assert that we have appropriate pigment if we have more than 3 cone types
            if (numel(obj.coneDensities)>3) && (any(obj.coneDensities(4:end)>0.0))
                assert(numel(obj.coneDensities) == size(obj.pigment.absorptance,3), ...
                    sprintf('cPhotoPigment is not initialized for %d types of cones', numel(obj.coneDensities)));
            end
            
            % These listeners make sure the wavelength support
            % in obj.pigment and obj.macular match the wave property
            addlistener(obj.pigment, 'wave', 'PostSet', @obj.matchWaveInAttachedMacular);
            addlistener(obj.macular, 'wave', 'PostSet', @obj.matchWaveInAttachedPhotopigment);
            
            
            % These listeners trigger the assignConeTypes() method when the
            % coneDensities or the tritanopicRadiusDegs properties are
            % assigned by the user.
            addlistener(obj, 'coneDensities','PostSet', @obj.assignConeTypes);
            addlistener(obj, 'tritanopicRadiusDegs', 'PostSet', @obj.assignConeTypes);
                
            if (isempty(p.Results.coneData))
                if (p.Results.computeMeshFromScratch)
                    % Re-generate lattice
                    obj.regenerateConePositions(...
                        p.Results.maxMeshIterations,  ...
                        p.Results.visualizeMeshConvergence, ...
                        p.Results.exportMeshConvergenceHistoryToFile, ...
                        'customMinRFspacing', customMinRFspacing, ...
                        'customDegsToMMsConversionFunction', customDegsToMMsConversionFunction, ...
                        'customMMsToDegsConversionFunction', customMMsToDegsConversionFunction, ...
                        'customRFspacingFunction', customRFspacingFunction);
                else
                    % Import positions by cropping a large pre-computed patch
                    obj.initializeConePositions();
                end

                % Remove cones within the optic disk
                obj.removeConesWithinOpticNerveHead();

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
            

            % Compute min and max cone position
            obj.minRFpositionMicrons = squeeze(min(obj.coneRFpositionsMicrons,[],1));
            obj.maxRFpositionMicrons = squeeze(max(obj.coneRFpositionsMicrons,[],1));
            
            % Compute ecc in microns
            if (isempty(obj.minRFpositionMicrons))
                % No cones, set value differently
                if (~isempty(obj.micronsPerDegreeApproximation))
                    obj.eccentricityMicrons = obj.eccentricityDegs * obj.micronsPerDegreeApproximation;
                else
                    obj.eccentricityMicrons = obj.eccentricityDegs *  300;
                end
            else
                obj.eccentricityMicrons = 0.5*(obj.minRFpositionMicrons + obj.maxRFpositionMicrons);
            end
            
    
            % Compute photon absorption attenuation factors to account for
            % the decrease in outer segment legth with ecc.
            obj.computeOuterSegmentLengthEccVariationAttenuationFactors('useParfor', obj.useParfor);
        end
        
        % Method to visualize the cone mosaic and its activation
        visualize(obj, varargin);
        
        % Method to visualize the continuous, full absorptions density and the actual cone positions
        visualizeFullAbsorptionsDensity(obj, figNo);
        
        % Method to generate eye movement sequences
        emGenSequence(obj, durationSeconds, varargin);
        
        % Method to generate an ensemble of OIs for the mosaic
        [oiEnsemble, psfEnsemble, zCoeffs] = oiEnsembleGenerate(obj, oiSamplingGridDegs, varargin);
        
        % Method to compute the mosaic response
        [absorptionsCount, noisyAbsorptionInstances, ...
            photoCurrents, photoCurrentInstances, responseTemporalSupport]  = compute(obj, oi, varargin);
        
        % Method for suggesting a minimal pixel size (in degrees) for the input scene
        % depending on the aperture size of the mosaic's cones and whether
        % the computation is done using ecc-dependent blur mode
        scenePixelSizeDegs = suggestedScenePixelSizeDegs(obj, eccVaryingConeBlur);
        
        % Method to return indices of cones within an ROI
        coneIndices = indicesOfConesWithinROI(obj, roi);
        
        % Generate struct representing the optical disk
        [odStructMicrons, odStructDegs] = odStruct(obj);
        
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
            if (isempty(obj.coneRFpositionsMicrons))
                % No cones in the current patch. Select some value
                if (~isempty(obj.micronsPerDegreeApproximation))
                    val = obj.micronsPerDegreeApproximation;
                else
                    val = 300;
                end
            else
                val = (max(obj.coneRFpositionsMicrons,[],1) - min(obj.coneRFpositionsMicrons,[],1)) / ...
                      (max(obj.coneRFpositionsDegs,[],1)    - min(obj.coneRFpositionsDegs,[],1));
            end
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
        regenerateConePositions(obj, maxIterations, visualizeConvergence, exportHistoryToFile, varargin);
        
        % Remove cones located within the optic disk
        removeConesWithinOpticNerveHead(obj);
        
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

        % Called automatically after either the coneDensity of the tritanopicRadiusDegs are set
        assignConeTypes(obj, src, ~);
        
        
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
        
        % Static method to validate an ROI struct
        validateROI(roi);
    
        % Static method to generate an ROI outline from an ROIstruct
        roiOutline = generateOutline(roi);
    
        % Static method to convert an ROIoutline in degs to an ROIoutline in microns
        roiOutlineMicrons = convertOutlineToMicrons(roiOutlineDegs,micronsPerDegreeApproximation);
        
        % Compute a 2D cone density map
        coneDensityMap = densityMap(rfPositions,rfSpacings, sampledPositions);
        
        % Function to generate a semitransparent controur plot
        semiTransparentContourPlot(axesHandle, xSupport, ySupport, zData, zLevels, cmap, alpha, contourLineColor);
    end
end

