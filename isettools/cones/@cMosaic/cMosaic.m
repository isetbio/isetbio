classdef cMosaic < handle
    % Create an optimized cone mosaic object
    %
    % Syntax:
    %   theConeMosaic = cMosaic;
    %   theConeMosaic = cMosaic('pigment', pp);
    %   theConeMosaic = cMosaic('sizeDegs',[5,5],'eccentricityDegs',[0,0]);
    %
    % Description:
    %    An object of the @cMosaic class contains parameters to create a
    %    retinal cone mosaic and enable computation of cone excitations and
    %    photocurrent response.
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
    %    'rodIntrusionAdjustedConeAperture' - Logical, indicating whether to adjust the cone aperture (light collecting area), or scalar in (0..1], to account for rods. Default: false 
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
        
        % Boolean indicating whether to adjust cone aperture to account for
        % rod intrusion - in effect only if eccVaryingConeAperture is true
        rodIntrusionAdjustedConeAperture;

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
        
        % No ecc-varying correction for any parameter
        anchorAllEccVaryingParamsToTheirFovealValues;

        % Integration time (seconds)
        integrationTime;
        
        % Cone diameter to cone spacing ratio (1.0 for compatibility with
        % old cone mosaics). But this should be < 1.0 to allow
        % for rod spacing. According to Chen et al (1993), Figure 4, cone
        % spacing is 2.8 microns and cone inner segment diameter is 2.3, so an
        % appropriate value for this coneDiameterToSpacingRatio would be
        % 2.3/2.8 = 0.821
        coneDiameterToSpacingRatio = 1.0;
        
        % Custom cone aperture params (struct)
        coneApertureModifiers = [];
        
        % Factor for cone coupling (the spatial constant of cone coupling
        % will equal this factor x cone diameter
        coneCouplingLambda;
        
        % Poisson noise flag for cone excitations 
        noiseFlag;

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
        
        % coupling weights - computed the first time compute() is called
        coneCouplingWeights;
        neighboringCoupledConeIndices;
        
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
        
        % Position [x,y] of the center of the cone mosaic (degs)
        % The variable should be renamed positionDegs, but ...
        eccentricityDegs;
        
        % Eccentricity [x,y] of the center of the cone mosaic (microns)
        % Should be positionMicrons
        eccentricityMicrons; 
        
        % Retinal quadrant (i.e. nasal, temporal, inferior, superior
        horizontalRetinalMeridian;
        verticalRetinalMeridian;

        % Size of the cone mosaic [width, height]
        sizeDegs;
        
        % Either 'right eye' or 'left eye'
        whichEye;
        
        % Fixational eye movement object for the mosaic
        fixEMobj = [];
        
        % Photon absorption attenuation factors to account for decrease in 
        % outer segment length with eccentricity
        outerSegmentLengthEccVariationAttenuationFactors = [];

        % Cone aperture diameter (summation) for each cone
        coneApertureDiametersMicrons = [];
        coneApertureDiametersDegs = [];
        
        % This depends on the value of the apertureModifier
        coneApertureToDiameterRatio = [];
        
        % This depends on the value of the apertureModifier
        coneApertureToConeCharacteristicRadiusConversionFactor = [];

        % Resolution with which zoning is performed
        oiResMicronsForZoning = [];
        
        % Indices of cones in each zone
        coneIndicesInZones = [];
        
        % Blur diameter for the different zones of cones
        blurApertureDiameterMicronsZones = [];
        blurApertureDiameterDegsZones = []; 
        
        % User-settable microns per degree
        micronsPerDegreeApproximation = [];
        
        % Custom degs <-> mm conversions
        customDegsToMMsConversionFunction = [];
        customMMsToDegsConversionFunction = [];
            
        % Achieved cone densities
        achievedConeDensities;

        % Flag indicating whether the cMosaic is based on i mported conedata
        employsImportedConeData = false;
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

        % Number of cones
        conesNum;
    end
    
    
    % Private properties
    properties (GetAccess=private, SetAccess=private)
        % Size (in degs) of source lattice from which to crop positions for
        % the desired eccentricity
        sourceLatticeSizeDegs = 58;
        
        % Scalar or empty, indicating whether to check and eliminate cones that are
        % too close to each other (overlap = overlappingConeFractionForElimination). 
        % Cone overlap is rare, but it happens at
        % some positions within the pre-computed mosaic. If this flag is
        % enabled during instantiation, we search for overlapping cones and
        % eliminate them. Depending on the size of the mosaic, this process
        % may take a while.
        overlappingConeFractionForElimination = [];

        % Struct containing information related to the current cone partition
        % into zones based on their aperture size
        cachedConePartition = [];
        
        % Struct containing information related to the current
        % macularPigmentDensityBoostFactors
        cachedMacularPigmentDensityBoostFactors = [];
        
        % OSLength attenuation factors (only used when importing coneData)
        importedOSLengthAttenuationFactors = [];
        
        % Cone aperture shrinkage factors due to progressive rod intrusion 
        % with eccentricity
        coneApertureRodIntrusionInducedShrinkageFactors = [];

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
        
        % The theoretical cone densities
        theoreticalConeDensities;
    end
    
    % Public methods
    methods
        
        % Constructor
        function [obj, cmParams] = cMosaic(varargin)
            
            cmParams = cMosaicParams;  % Default.  Updated at the end
            
            % Parse input
            p = inputParser;
            p.addParameter('name', 'cone mosaic', @ischar);
            p.addParameter('sourceLatticeSizeDegs', 58, @isscalar);
            p.addParameter('overlappingConeFractionForElimination', [], @(x)( (isempty(x)) || (isscalar(x)&&(x<=1.0)&&(x>0))));
            p.addParameter('wave', 400:10:700, @isnumeric);
            p.addParameter('pigment', cPhotoPigment(), @(x) isa(x, 'cPhotoPigment'));
            p.addParameter('macular', Macular(), @(x)isa(x, 'Macular'));
            p.addParameter('coneData', [], @(x)(isempty(x) || (isstruct(x))));
            
            % These are synonyms.  If positionDegs is sent in, it
            % overwrites eccentricityDegs
            p.addParameter('eccentricityDegs', [], @(x)(isnumeric(x) && (numel(x) == 2)));
            p.addParameter('positionDegs', [], @(x)(isnumeric(x) && (numel(x) == 2)));
            
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
            p.addParameter('rodIntrusionAdjustedConeAperture', false, @(x) ((islogical(x))||((isscalar(x))&&((x>0)&&(x<=1)))));
            p.addParameter('eccVaryingConeAperture', true, @islogical);
            p.addParameter('eccVaryingConeBlur', false, @islogical);
            p.addParameter('eccVaryingOuterSegmentLength', true, @islogical);
            p.addParameter('eccVaryingMacularPigmentDensity', true, @islogical);
            p.addParameter('eccVaryingMacularPigmentDensityDynamic', false, @islogical);
            p.addParameter('anchorAllEccVaryingParamsToTheirFovealValues', false, @islogical);
            p.addParameter('coneCouplingLambda', [], @(x)(isempty(x)||isscalar(x)));
            p.addParameter('coneApertureModifiers', struct('smoothLocalVariations', true), @isstruct);
            p.addParameter('coneDiameterToSpacingRatio', 1.0,  @(x)(isscalar(x)&&(x<=1.0)));
            p.addParameter('coneDensities', [0.6 0.3 0.1 0.0], @(x)(isnumeric(x) && ((numel(x) == 3)||(numel(x)==4))));
            p.addParameter('tritanopicRadiusDegs', 0.15, @isscalar);
            p.addParameter('integrationTime', 5/1000, @isscalar);
            p.addParameter('opticalImagePositionDegs', 'mosaic-centered', @(x)(ischar(x) || (isnumeric(x)&&numel(x)==2)));
            p.addParameter('useParfor', true, @islogical);
            
            % Can we have frozen noise without a randomseed? And, if we
            % have a seed, perhaps we intend frozen noise? So we agree on
            % what we want and deal with it here and then below around line
            % 490.
            p.addParameter('noiseFlag', 'random', @(x)(ischar(x) && (ismember(x, {'random', 'frozen', 'none'}))));
            p.addParameter('randomSeed', [], @(x)(isempty(x) || isscalar(x)));

            p.parse(varargin{:});
            
            obj.name    = p.Results.name;
            obj.sourceLatticeSizeDegs = p.Results.sourceLatticeSizeDegs;
            obj.overlappingConeFractionForElimination = p.Results.overlappingConeFractionForElimination;
            obj.macular = p.Results.macular;
            obj.pigment = p.Results.pigment;
            obj.wave    = p.Results.wave;
            
            % BW wants to use positionDegs rather than eccentricity.  For
            % teaching purposes.
            if ~isempty(p.Results.eccentricityDegs)
                % eccentricityDegs is always used if present, for
                % historical reasons 
                obj.eccentricityDegs = p.Results.eccentricityDegs;
            elseif ~isempty(p.Results.positionDegs)
                % Use positionDegs if available but there is not
                % eccentricityDegs
                obj.eccentricityDegs = p.Results.positionDegs;
            else                
                % Default
                obj.eccentricityDegs = [0,0];
            end
                
            obj.sizeDegs = p.Results.sizeDegs;
            obj.whichEye = p.Results.whichEye;
            
            obj.rodIntrusionAdjustedConeAperture = p.Results.rodIntrusionAdjustedConeAperture;
            obj.eccVaryingConeAperture = p.Results.eccVaryingConeAperture;
            obj.eccVaryingOuterSegmentLength = p.Results.eccVaryingOuterSegmentLength;
            obj.eccVaryingMacularPigmentDensity = p.Results.eccVaryingMacularPigmentDensity;
            obj.eccVaryingConeBlur = p.Results.eccVaryingConeBlur;
            obj.eccVaryingMacularPigmentDensityDynamic  = p.Results.eccVaryingMacularPigmentDensityDynamic;
            obj.anchorAllEccVaryingParamsToTheirFovealValues = p.Results.anchorAllEccVaryingParamsToTheirFovealValues;

            obj.coneCouplingLambda = p.Results.coneCouplingLambda;
            
            obj.coneDensities = p.Results.coneDensities;
            obj.tritanopicRadiusDegs = p.Results.tritanopicRadiusDegs;
            obj.coneApertureModifiers = p.Results.coneApertureModifiers;
            obj.coneDiameterToSpacingRatio = p.Results.coneDiameterToSpacingRatio;
            
            obj.noiseFlag  = p.Results.noiseFlag;
            obj.randomSeed = p.Results.randomSeed;
            obj.integrationTime = p.Results.integrationTime;
            obj.opticalImagePositionDegs = p.Results.opticalImagePositionDegs;
            
            % Parallel computations
            obj.useParfor = p.Results.useParfor;
            
            % Custom mm/deg conversion factors
            obj.customDegsToMMsConversionFunction = p.Results.customDegsToMMsConversionFunction;
            obj.customMMsToDegsConversionFunction = p.Results.customMMsToDegsConversionFunction;
            obj.micronsPerDegreeApproximation = p.Results.micronsPerDegree;
            
            % Assert that only one microns/deg conversion has been specified
            if (...
               (~isempty(obj.micronsPerDegreeApproximation)) && ...
               ((~isempty(obj.customDegsToMMsConversionFunction)) || (~isempty(obj.customMMsToDegsConversionFunction)))...
               )
               error('A custom ''micronsPerDegree'' and a deg <-> MM conversion function cannot be both specified.');
            end
            
            % Assert that if one conversion function is specified its
            % reverse is also specified
            if ((isempty(obj.customDegsToMMsConversionFunction)) && (~isempty(obj.customMMsToDegsConversionFunction)))
               error('A ''customMMsToDegsConversionFunction'' was specified but not its inverse');
            end
            
            if ((~isempty(obj.customDegsToMMsConversionFunction)) && (isempty(obj.customMMsToDegsConversionFunction)))
               error('A ''customDegsToMMsConversionFunction'' was specified but not its inverse');
            end
            
            % Assert that we have appropriate pigment if we have more than 3 cone types
            if (numel(obj.coneDensities)>3) && (any(obj.coneDensities(4:end)>0.0))
                assert(numel(obj.coneDensities) == size(obj.pigment.absorptance,3), ...
                    sprintf('cPhotoPigment is not initialized for %d types of cones', numel(obj.coneDensities)));
            end

            % Assert that if anchorAllEccVaryingParamsToTheirFovealValues
            % is set, all ecc-dependent flag are off
            if (obj.anchorAllEccVaryingParamsToTheirFovealValues)
                assert(obj.eccVaryingConeAperture == false, ...
                    sprintf('Can not have both ''anchorAllEccVaryingParamsToTheirFovealValues'' and ''eccVaryingConeAperture'' set to true.'));
                assert(obj.eccVaryingOuterSegmentLength == false, ...
                    sprintf('Can not have both ''anchorAllEccVaryingParamsToTheirFovealValues'' and ''eccVaryingOuterSegmentLength'' set to true.'));
                assert(obj.eccVaryingMacularPigmentDensity == false, ...
                    sprintf('Can not have both ''anchorAllEccVaryingParamsToTheirFovealValues'' and ''eccVaryingMacularPigmentDensity'' set to true.'));
                assert(obj.eccVaryingMacularPigmentDensityDynamic == false, ...
                    sprintf('Can not have both ''anchorAllEccVaryingParamsToTheirFovealValues'' and ''eccVaryingMacularPigmentDensityDynamic'' set to true.'));
                assert(obj.eccVaryingConeBlur == false, ...
                    sprintf('Can not have both ''anchorAllEccVaryingParamsToTheirFovealValues'' and ''eccVaryingConeBlur'' set to true.'));
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
                    
                    % Custom mesh generation function 
                    customMinRFspacing = p.Results.customMinRFspacing;
                    customRFspacingFunction = p.Results.customRFspacingFunction;
            
                    % Re-generate lattice
                    obj.regenerateConePositions(...
                        p.Results.maxMeshIterations,  ...
                        p.Results.visualizeMeshConvergence, ...
                        p.Results.exportMeshConvergenceHistoryToFile, ...
                        'customMinRFspacing', customMinRFspacing, ...
                        'customRFspacingFunction', customRFspacingFunction);
                else
                    % Import positions by cropping a large pre-computed patch
                    obj.initializeConePositions(obj.overlappingConeFractionForElimination);
                end

                % Remove cones within the optic disk
                obj.removeConesWithinOpticNerveHead();

                % Set random seed so we get consistent cone assignment
                if (~isempty(obj.randomSeed))
                    rng(obj.randomSeed);
                end

                % Assign cone types
                obj.assignConeTypes();

                % If we want random responses, but we were not passed a
                % random seed, make a noise seed
                if isequal(obj.noiseFlag,'random') && (isempty(obj.randomSeed))
                    % Make a noise seed and store it.  One in a million.
                    obj.randomSeed = randi(1e6);

                elseif isequal(obj.noiseFlag,'frozen') && isempty(obj.randomSeed)
                    error('Frozen noise but no seed provided.');
                end

                
            else
                % Load cone positions and cone types from passed struct
                obj.importExternalConeData(p.Results.coneData);
                
                % Set random seed
                if (isempty(obj.randomSeed))
                    rng('shuffle');   % Initialize rng with current time
                    
                    % Remember the seed in case we freeze the noise next
                    % round.
                    obj.randomSeed = randi(1e4);
                else
                    disp(obj.randomSeed)
                end
                if (~isempty(obj.randomSeed))
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
                elseif (~isempty(obj.customDegsToMMsConversionFunction))
                    obj.eccentricityMicrons = obj.customDegsToMMsConversionFunction(obj.eccentricityDegs)*1e3;
                else
                    error('The cMosaic has no cones, no specification for ''micronsPerDegreeApproximation'' and no specification for ''obj.customDegsToMMsConversionFunction''. Unable to determine eccentricity in microns.')
                end

            else
                obj.eccentricityMicrons = 0.5*(obj.minRFpositionMicrons + obj.maxRFpositionMicrons);
            end
            
            % Set the corresponding horizontal and vertical retinal meridian 
            [obj.horizontalRetinalMeridian, obj.verticalRetinalMeridian] = cMosaic.retinalMeridiansForEccentricityInEye(...
                obj.eccentricityDegs(1), obj.eccentricityDegs(2),  obj.whichEye);


            % Compute photon absorption attenuation factors to account for
            % the decrease in outer segment legth with ecc.
            obj.computeOuterSegmentLengthEccVariationAttenuationFactors('useParfor', obj.useParfor);
            
            
            if nargout > 1
                % Return a struct with all the Results parameters from this
                % object. If you want the default, use
                %
                %   cmParams = cMosaicParams;
                %
                % But that needs people to update if they change the
                % parameters in this routine.
                cmParams = p.Results;
                
                % Sometimes people send in positionDegs instead of
                % eccentricityDegs (BW).  But eccentricityDeg cannot be
                % empty.  So we fill it here.
                cmParams.eccentricityDegs = obj.eccentricityDegs;
                
                % A random seed was used.  We store it in the return in
                % case we want it in the future
                cmParams.randomSeed = obj.randomSeed;
            end
            
        end % constructor
        
        % Method to crop the cMosaic to only incluce cones with specific
        % indices
        cropMosaicToIncluceConesWithIndices(obj,keptConeIndices);
        
        % Method to transform a distance specified in units of retinal
        % microns to a distance specified in units of visual degrees based on the @cMosaic configuration
        microns = distanceDegreesToDistanceMicronsForCmosaic(obj, degrees);
        
        % Method to transform a distance specified in units of visual
        % degrees to a distance specified in units of retinal microns based on the @cMosaic configuration
        degrees = distanceMicronsToDistanceDegreesForCmosaic(obj, microns);
        
        % Method to transform a SIZE specified in units of retinal microns at a
        % given eccentricity (also specified in microns) to a SIZE in units of 
        % visual degrees based on the @cMosaic configuration
        sizeDegrees = sizeMicronsToSizeDegreesForCmosaic(obj, sizeMicrons, eccentricityMicrons);

        % Method to transform a SIZE specified in units of visual degrees at a
        % given eccentricity (also specified in degrees) to a SIZE in units of 
        % retinal microns based on the @cMosaic configuration
        sizeMicrons = sizeDegreesToSizeMicronsForCmosaic(obj, sizeDegrees, eccentricityDegrees);
    
        % Method for generating the cone aperture blur kernel
        apertureKernel = generateApertureKernel(obj, coneApertureDiameterMicrons, oiResMicrons, lowOpticalImageResolutionWarning);
        
        % Blur sigma (in microns) of a cone with index theConeIndex, from its blur zone
        % Only for Gaussian aperture modifier
        apertureBlurSigmaMicrons = apertureBlurSigmaMicronsOfConeFromItsBlurZone(obj, theConeIndex);
        
        % Method to compute aperture areas
        apertureAreasMetersSquared = computeApertureAreasMetersSquared(obj, coneIndices);

        % Method to compute the effective os length attenuation factors
        v = computeEffectiveOSlengthAttenuationFactors(obj, coneIndices);

        % Method to visualize the cone mosaic and its activation
        params = visualize(obj, varargin);
        
        % Method to visualize activation profiles for particular cone types along a
        % horizontal slice ROI
        hFig = visualizeHorizontalConeActivationProfiles(obj, theConeMosaicResponse, coneTypesToVisualize, ...
            horizontalSliceYcoordDegs, horizontalSliceWidthDegs, maxResponse, varargin);



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
        
        % Method to return indices of cones within a geometry struct appropriate for @regionOfInterest
        coneIndices = indicesOfConesWithinROI(obj, geometryStruct);
        
        % Method to return the cone indices of cones of a certain type that
        % lie within a horizontal rect ROI
        [coneIndicesOfCertainConeType, theROI] = ...
            coneIndicesOfCertainTypeWithinHorizontalSlice(obj, slicePositionDegs, sliceWidthDegs, coneType);

        % Method to reassigne the type of any set of cones, specified by
        % their cone index
        reassignTypeOfCones(obj, coneIndices, newConeType);

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

        function val = get.conesNum(obj)
            val = size(obj.coneRFpositionsMicrons,1);
        end

        function set.coneDensities(obj, val)
            obj.theoreticalConeDensities = val;
        end
        
        function val = get.coneDensities(obj)
            if (isempty(obj.lConeIndices))
                val = obj.theoreticalConeDensities;
            else
                % Return actual densities
                conesNum = size(obj.coneRFpositionsMicrons,1);
                val = [...
                    numel(obj.lConeIndices)/conesNum ...
                    numel(obj.mConeIndices)/conesNum ...
                    numel(obj.sConeIndices)/conesNum ...
                    numel(obj.kConeIndices)/conesNum];
            end
        end
        
        function set.eccVaryingConeAperture(obj, val)
            obj.eccVaryingConeAperture = val;
            if (~isempty(obj.coneRFpositionsMicrons))
                lowOpticalImageResolutionWarning = true;
                obj.computeConeApertures(lowOpticalImageResolutionWarning);
            end
        end
        
        function set.rodIntrusionAdjustedConeAperture(obj,val)
            obj.rodIntrusionAdjustedConeAperture = val;
            if (~isempty(obj.coneRFpositionsMicrons))
                lowOpticalImageResolutionWarning = true;
                obj.computeConeApertures(lowOpticalImageResolutionWarning);
            end
        end

        function set.eccVaryingConeBlur(obj, val)
            obj.eccVaryingConeBlur = val;
            if (~isempty(obj.coneRFpositionsMicrons))
                lowOpticalImageResolutionWarning = true;
                obj.computeConeApertures(lowOpticalImageResolutionWarning);
            end
        end
        
        
        function set.coneApertureModifiers(obj, val)
            % Validate fields
            fNames = fieldnames(val);
            validFieldNames = {'smoothLocalVariations', 'shape', 'sigma', 'PillboxApertureIntegrationMatching'};
            for k = 1:numel(fNames)
                assert(ismember(fNames{k}, validFieldNames), ...
                    sprintf('Unknown field ''%s'' in coneApertureModifier struct', fNames{k}));
            end
            
            % Cone aperture (light collecting disk) has a diameter that 
            % is 0.79 x cone diameter. We chose this to be consistent with
            % coneMosaicHex, .e.g.: 
            % c = coneMosaicHex(3);
            % c.pigment.pdWidth/c.pigment.width, which gives 0.7899.
            %
            % Note. The coneApertureToDiameterRatio is overriden when a
            % coneApertureModifier struct is passed (see below) with a
            % a 'shape' field that is set to 'Gaussian'. In that case
            % it becomes equal to 1 (cMosaic.computeConeApertures.m, line 30)
            % and the Gaussian sigma parameter is set to control the aperture
            pillboxConeApertureToDiameterRatioToMatchConeMosaixHexAperture = 0.79;
        
            % If we want the pillbox aperture to have an integral that matches 
            % the integral of Gaussian area with sigma = 0.204 * aperture,
            % (as in 'Serial Spatial Filters in Vision', by Chen, Makous and Williams, 1993) 
            % set this to 2*sqrt(2)*0.204 = 0.577
            pillboxConeApertureToDiameterRatioToMatchCMW1993Aperture = 2*sqrt(2)*0.204;
            
            % If we have a Pillbox aperture shape, or if there is no
            % aperture shape parameter, determine the value of obj.coneApertureToDiameterRatio
            if (isfield(val, 'shape')) && (strcmp(val.shape,'Pillbox')) && (isfield(val, 'PillboxApertureIntegrationMatching'))
                switch (val.PillboxApertureIntegrationMatching)
                    case 'ChenMakousWilliams1993'
                        obj.coneApertureToDiameterRatio = pillboxConeApertureToDiameterRatioToMatchCMW1993Aperture;
                    case 'ConeMosaicHex'
                        obj.coneApertureToDiameterRatio = pillboxConeApertureToDiameterRatioToMatchConeMosaixHexAperture;
                end
            else
                obj.coneApertureToDiameterRatio = pillboxConeApertureToDiameterRatioToMatchConeMosaixHexAperture;
            end
            
            obj.coneApertureModifiers = val;
            if (~isempty(obj.coneRFpositionsMicrons))
                lowOpticalImageResolutionWarning = true;
                obj.computeConeApertures(lowOpticalImageResolutionWarning);
            end
        end
        
        function set.coneDiameterToSpacingRatio(obj,val)
            if (val <= 1.0)
                obj.coneDiameterToSpacingRatio = val;
                if (~isempty(obj.coneRFpositionsMicrons))
                    lowOpticalImageResolutionWarning = true;
                    obj.computeConeApertures(lowOpticalImageResolutionWarning);
                end
            else
                error('coneDiameterToSpacingRatio must be <= 1.0.');
            end
        end
        
        function set.coneCouplingLambda(obj, val)
            obj.coneCouplingLambda = val;
            obj.coneCouplingWeights = [];
            obj.neighboringCoupledConeIndices = [];
        end                

    end % Public methods
    
    
    methods (Access=private)
        % Initialize cone positions by importing them from a large previously-computed mesh
        initializeConePositions(obj, overlappingConeFractionForElimination);

        % Initialize cone positions by regenerating a new mesh. Can be slow.
        regenerateConePositions(obj, maxIterations, visualizeConvergence, exportHistoryToFile, varargin);
        
        % Remove cones located within the optic disk
        removeConesWithinOpticNerveHead(obj);
        
        % Crop data for desired ROI
        cropMosaicDataForDesiredROI(obj);
        
        % Update state of cMosaic given kept cone indices
        updateStateGivenKeptConeIndices(obj, keptConeIndices);
        
        % Function to re-compute the cone apertures
        computeConeApertures(obj, lowOpticalImageResolutionWarning);
     
        % Method to electrically couple cone responses
        coupledResponses = electricallyCoupleConeResponses(obj, responses);
        
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
        % Static method to compute the macular pigment density boost
        % factors at an arbitray array of retinal positions (in degrees)
        macularPigmentDensityBoostFactors = macularPigmentBoostFactors(theMacular, retinalPositionsDegs)

        % Return the outer segment length in microns at the passed
        % eccentricity as well as the foveal outer segment length using
        % data from Banks, Sekuler and Anderson (1991). "Peripheral spatial vision: limits
        % imposed by optics, photoreceptors and receptor pooling".
        [osLengthMicrons, osLengthMicronsFoveal] = outerSegmentLengthFromEccentricity(eccDegs);

        % Static method to return signed horizontal and vertical
        % eccentricities corresponding to radial eccentricities specified
        % on one of the 4 principal retinal meridians and eye
        [horizontalEcc, verticalEcc] = eccentricitiesForRetinaMeridianInEye(...
            radialEcc, retinaQuadrant, whichEye);

        % Static method to return retinal meridians corresponding to
        % signed horizontal and vertical eccentricities in a specified eye
        [horizontalRetinalMeridian, verticalRetinalMeridian] = ...
            retinalMeridiansForEccentricityInEye(horizontalEcc, verticalEcc, whichEye);

        % Static method to generate noisy absorption response instances
        % from the mean absorption responses
        noisyAbsorptionInstances = noisyInstances(meanAbsorptions, varargin);
        
        % Compute a 2D cone density map
        coneDensityMap = densityMap(rfPositions,rfSpacings, sampledPositions);
        
        % Function to generate a semitransparent controur plot
        semiTransparentContourPlot(axesHandle, xSupport, ySupport, zData, zLevels, cmap, alpha, contourLineColor, varargin);
    
        % Function for identifying overlapping RFs
        [rfsToKeep, rfsToBeEliminated, overlappingOtherRFs] = identifyOverlappingRFs(xPos, yPos, ...
             RFpositionsMicrons, RFspacingsMicrons, maxSeparationForDeclaringOverlap);
    end
end

