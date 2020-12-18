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
    %    12/../2020  NPC  Wrote it
    %    12/18/20  dhb    Use new quantalEfficiency property of photoPigment.

    
    properties (Constant)
        % Cone types
        LCONE_ID = 1;
        MCONE_ID = 2;
        SCONE_ID = 3;
        KCONE_ID = 4;
        
        % Assume cone aperture (light collecting disk) has a diameter that 
        % is 0.7 x cone spacing
        coneApertureToDiameterRatio = 0.7;
    end
    
    % Public properties
    properties  (GetAccess=public, SetAccess=public)
        % Name of the mosaic
        name;
        
        % Cone photopigment for the mosaic
        pigment;
        
        % Macular pigment object for mosaic
        macular;
        
        % Eccentricity [x,y] of the center of the cone mosaic 
        eccentricityDegs;
        
        % Size of the cone mosaic [width, height]
        sizeDegs;
        
        % Either 'right eye' or 'left eye'
        whichEye;

        % Whether blurring by cone aperture will be performed using the
        % median cone aperture (false, less accurate, faster), or using 
        % the eccentricity-varying cone aperture (true, more accurate, slower)
        eccVaryingConeBlur;

        % 3- or 4-element vector containing relative cone densities
        coneDensities;
        
        % Radius of S-cone free area around (0,0)
        tritanopicRadiusDegs;
        
        % Poisson noise flag for cone excitations 
        noiseFlag;
        
        % Color for cone rendering in mosaic visualization
        lConeColor = [1 0.5 0.6];
        mConeColor = [0.2 1 0.5];
        sConeColor = [0.5 0.1 1];
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
    end
    
    
    % Private properties
    properties (GetAccess=private, SetAccess=private)
        % Size (in degs) of source lattice from which to crop positions for
        % the desired eccentricity
        sourceLatticeSizeDegs = 45;

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
            p.addParameter('eccentricityDegs', [0 0], @(x)(isnumeric(x) && (numel(x) == 2)));
            p.addParameter('sizeDegs', [0.2 0.2], @(x)(isnumeric(x) && (numel(x) == 2)));
            p.addParameter('whichEye', 'right eye', @(x)(ischar(x) && (ismember(x, {'left eye', 'right eye'}))));
            p.addParameter('eccVaryingConeBlur', false, @islogical);
            p.addParameter('coneDensities', [0.6 0.3 0.1 0.0], @(x)(isnumeric(x) && ((numel(x) == 3)||(numel(x)==4))));
            p.addParameter('tritanopicRadiusDegs', 0.15, @isscalar);
            p.addParameter('noiseFlag', 'random', @(x)(ischar(x) && (ismember(x, {'random', 'frozen', 'none'}))));
            p.parse(varargin{:});
            
            obj.name = p.Results.name;
            obj.macular = p.Results.macular;
            obj.pigment = p.Results.pigment;
            obj.wave = p.Results.wave;
            
            obj.eccentricityDegs = p.Results.eccentricityDegs;
            obj.sizeDegs = p.Results.sizeDegs;
            obj.whichEye = p.Results.whichEye;

            obj.eccVaryingConeBlur = p.Results.eccVaryingConeBlur;

            obj.coneDensities = p.Results.coneDensities;
            obj.tritanopicRadiusDegs = p.Results.tritanopicRadiusDegs;
            obj.noiseFlag = p.Results;
            
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
            
            % Initialize positions
            obj.initializeConePositions();
            
            % Assign types
            obj.assignConeTypes();
        end
        
        % Method to assign cone types
        assignConeTypes(obj, varargin);
        
        % Method to visualize the cone mosaic
        visualize(obj, varargin);
        
        % Getter/Setter methods for dependent variables
        % QE
        function val = get.qe(obj)
            % Compute the quantal efficiency of the cones, incorporating
            % macular pigment.      
            val = bsxfun(@times, obj.pigment.quantalEfficiency, ...
                obj.macular.transmittance) ;
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
        
    end
    
    methods (Access=private)
        initializeConePositions(obj);
        
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
        m = displayAbsorptionsDensity(ax, oi, absorptionsDensity, meanLevel, superimposeHorizontalSlice);
    end
    
end

