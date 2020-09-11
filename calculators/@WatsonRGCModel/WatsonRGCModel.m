classdef WatsonRGCModel
    % Create a WatsonRGCModel
    %
    %
    % Syntax:
    %   w = WatsonRGCModel('generateAllFigures', false);
    %
    % Usage:
    % - Compute peak cone density:
    %   peakConeDensityPerMM2 = WatsonRGCCalc.peakConeDensity('cones per mm2')   
    %   peakConeDensityPerMM2 = WatsonRGCCalc.peakConeDensity('cones per deg2')   
    %
    % References:
    %    Watson (2014). 'A formula for human RGC receptive field density as
    %    a function of visual field location', JOV (2014), 14(7), 1-17.
    %
    % History:
    %    11/8/19  NPC, ISETBIO Team     Wrote it.
    
    % Constant properties (model parameters)
    properties (Constant)
        % Cell array with meridian parameters
        % These meridians are in the Right Eye visual field domain See
        % Watson (2014) section titled "Conventions regarding meridians
        % ..." in the Introduction.
        
        % Watson (2014) orders the meridians temporal, superior, nasal, inferior, consistent with increasing polar angle 
        % in the visual field of the right eye
        
        enumeratedMeridianNames = {...
            'temporal meridian' ...
            'superior meridian' ...
            'nasal meridian' ...
            'inferior meridian' ...
            };
        
        enumeratedMeridianAngles = (0:3)*90;
        
        meridianParamsTable = {
            WatsonRGCModel.enumeratedMeridianNames{1}  struct('a_k', 0.9851, 'r_2k', 1.058,  'r_ek', 22.14); ... 
            WatsonRGCModel.enumeratedMeridianNames{2}  struct('a_k', 0.9935, 'r_2k', 1.035,  'r_ek', 16.35); ...
            WatsonRGCModel.enumeratedMeridianNames{3}  struct('a_k', 0.9729, 'r_2k', 1.084,  'r_ek',  7.633); ... 
            WatsonRGCModel.enumeratedMeridianNames{4}  struct('a_k', 0.996,  'r_2k', 0.9932, 'r_ek', 12.13) ... 
        };
    
        % Dictionary with meridian colors
        meridianColors = containers.Map(...
             { ...
             WatsonRGCModel.meridianParamsTable{1,1}, ...
             WatsonRGCModel.meridianParamsTable{2,1}, ...
             WatsonRGCModel.meridianParamsTable{3,1}, ...
             WatsonRGCModel.meridianParamsTable{4,1} ...
             }, ...
             { ...
             '[1.0 0.0 0.0]', ...
             '[0.0 0.0 1.0]', ...
             '[0.0 0.8 0.0]', ...
             '[0.2 0.2 0.2]', ...
             }, 'UniformValues', true);
            
        % Various acronyms and their meaning in the Watson (2014) paper
        glossaryTable = {
             'mRGCf'    'midget RGC receptive field'; ...
             'g'        'RGC'; ...
             'm'        'midget RGC'; ...
             'c'        'cone'; ...
             'gf'       'RGC receptive field'; ...
             'mf'       'midget RGC receptive field'; ...
             'dc(0)'    'peak cone density (at 0 deg eccentricity)'; ...
             'f0'       'fraction of all ganglion cells that are midget at 0 deg eccentricity'; ...
             'alpha',   'Conversion factor mm^2 -> deg^2 as a function of eccentricity' ...
        };
     
        % Fraction of all ganglion cells that are midget at 0 deg eccentricity
        f0 = 1/1.12;

        % Peak cone density (cones/deg^2) at 0 deg eccentricity (page 3, dc(0), Also in Appendix 4)
        dc0 = 14804.6;
        
        % Linear approximatio of microns/deg on the retina
        micronsPerDegreeLinearApproximation = 268;
            
        % Conversion factor, rho, of retinal distance deg->mm as as a function of eccentricity in
        % degs (Equation A5)
        rhoDegsToMMs = @(eccDegs) sign(eccDegs) .* (...
            0.268         * abs(eccDegs) + ...
            0.0003427     * abs(eccDegs) .^2 + ...
           -8.3309 * 1e-6 * abs(eccDegs) .^3);
        
        % Conversion factor, rho, of retinal distance mm->deg as as a function of eccentricity in
        % degs (Equation A6)
        rhoMMsToDegs = @(eccMM) sign(eccMM) .* (...
            3.556     * abs(eccMM) + ...
            0.05993   * abs(eccMM) .^2 + ...
           -0.007358  * abs(eccMM) .^3 + ...
            0.0003027 * abs(eccMM) .^4);
       
        % Conversion factor of size in visual degress at a given eccentricity in degrees to size in retinal microns 
        sizeDegsToSizeRetinalMicrons = @(sizeDegs, eccDegs) ...
            1000.0 * (WatsonRGCModel.rhoDegsToMMs(eccDegs + sizeDegs/2) - WatsonRGCModel.rhoDegsToMMs(eccDegs - sizeDegs/2))

        % Conversion factor of size in retinal microns at a given eccentricity in microns to size in visual degrees
        sizeRetinalMicronsToSizeDegs = @(sizeMicrons, eccMicrons) ...
            WatsonRGCModel.rhoMMsToDegs((eccMicrons + sizeMicrons/2)/1000.0) - WatsonRGCModel.rhoMMsToDegs((eccMicrons - sizeMicrons/2)/1000)
        
        % Conversion factor, alpha, of retinal area mm^2 -> deg^2 as a function of eccentricity in
        % degs (Equation A7)
        alpha = @(eccDegs) 0.0752 + ...
                    5.846 * 1e-5 * abs(eccDegs)    + ...
                   -1.064 * 1e-5 * abs(eccDegs).^2 + ...
                    4.116 * 1e-8 * abs(eccDegs).^3;
        
        % Convert density to spacing in a perfect hex mosaic. This is equation (A4) in the Watson (2014) paper.
        spacingFromDensity = @(density) sqrt(2.0./(sqrt(3.0)*density));

        % Convert density to spacing in a perfect hex mosaic. This is equation (A4) in the Watson (2014) paper.
        densityFromSpacing = @(spacing) 2.0./(sqrt(3.0)*spacing.^2);
        
        % Valid eccentricity units
        visualDegsEccUnits = 'deg';
        retinalMMEccUnits = 'mm';
        
        % Valid density units
        visualDegsDensityUnits = 'deg^2';
        retinalMMDensityUnits = 'mm^2';
        
        % Valid view names
        rightEyeVisualField = 'right eye visual field';
        rightEyeRetina = 'right eye retina';
        leftEyeRetina = 'left eye retina';
        
        % Ratio of cone aperture to cone-diameter, set to 0.7 in ISETBIO
        coneApertureToDiameterRatio = 0.7;
    end % Constant properties
    
    % Constant properties related to figure generation
    properties (Constant)
        paperTitleFull = 'Watson (2014): ''A formula for human RGC receptive field density as a function of visual field location'' ';
        paperTitleShort = 'Watson (2014) RGC model';
    end
   
    % Public properties (read-only)
    properties (SetAccess = private)
       % Dictionary with various acronyms of the the Watson (2014) paper and their meaning
       glossary;
        
       % Dictionary with meridian params indexed by meridian name
       meridianParams;
       
       % Struct with default preferences for all figures
       defaultFigurePrefs = struct(...
            'lineWidth', 1.5, ...
            'markerLineWidth', 1.0, ...
            'fontSize', 14, ...
            'fontAngle', 'italic', ...
            'grid', 'on', ...
            'backgroundColor', [1 1 1]);
    end
    
    % Public properties
    properties
        % Struct with default preferences for all figures
        figurePrefs
        
        % the default eccentricity support (in degs) for all figures
        eccDegs = 0:0.002:90;
    end
    
    % Public methods (class interface)
    methods
        % Constructor
        function obj = WatsonRGCModel(varargin) 
            % Parse input
            p = inputParser;
            p.addParameter('generateAllFigures', false, @islogical);
            p.addParameter('eccDegs', 0:0.002:90, @isnumeric);
            p.addParameter('instantiatePlotLab', false, @islogical);
            p.parse(varargin{:});
            
            % Set the default figure preferences
            obj.figurePrefs = obj.defaultFigurePrefs;
            
            % Set options
            generateAllFigures = p.Results.generateAllFigures;
            obj.eccDegs = p.Results.eccDegs;
            
            % Create dictionary with various acronyms of the the Watson (2014) paper and their meaning
            obj.glossary = containers.Map(obj.glossaryTable(:,1), obj.glossaryTable(:,2));
            
            % Create dictionary with meridian params 
            obj.meridianParams = containers.Map(obj.meridianParamsTable(:,1), obj.meridianParamsTable(:,2));

            if (p.Results.instantiatePlotLab)
                WatsonRGCModel.setUpPlotLab();
            end
            
            % Generate figures
            if (generateAllFigures)
                plotlabOBJ = WatsonRGCModel.unitTestFigure1('plotlabOBJ', []);
                plotlabOBJ = WatsonRGCModel.unitTestFigure5('plotlabOBJ', plotlabOBJ);
                plotlabOBJ = WatsonRGCModel.unitTestFigure9('plotlabOBJ', plotlabOBJ);
                plotlabOBJ = WatsonRGCModel.unitTestFigure10('plotlabOBJ', plotlabOBJ);
                plotlabOBJ = WatsonRGCModel.unitTestFigure11('plotlabOBJ', plotlabOBJ);
                plotlabOBJ = WatsonRGCModel.unitTestFigure14('plotlabOBJ', plotlabOBJ);
            end
        end % Constructor
        
        % Return cone RF spacing, cone RF density along the requested meridian and requested eccentricities
        [coneRFSpacing, coneRFDensity, rightEyeRetinalMeridianName] = coneRFSpacingAndDensityAlongMeridian(obj, eccentricities, rightEyeVisualFieldMeridianName, eccUnits, densityUnits, varargin);
        
        % Return midget RGC RF spacing and density along the requested meridian and requested eccentricities
        [mRGCRFSpacing, mRGCRFDensity, rightEyeRetinalMeridianName] = mRGCRFSpacingAndDensityAlongMeridian(obj, eccentricities, rightEyeVisualFieldMeridianName, eccUnits, densityUnits, varargin);
        
        % Return midget RFC RF spacing and density at each of the passed
        % retinal positions and the left/right eye
        [mRGCSpacings, mRGCRFDensities] = mRGCRFSpacingAndDensityAtRetinalPositions(obj, rfPositions, whichEye, posUnits, densityUnits, varargin);

        % Return cone RF spacing and density at each of the passed
        % retinal positions and the left/right eye
        [coneSpacings, coneDensities] = coneRFSpacingAndDensityAtRetinalPositions(obj, rfPositions, whichEye, posUnits, densityUnits, varargin);

        
        % Return total RGC RF density along the requested meridian and requested eccentricities
        totalRGCRFSpacing = totalRGCRFDensityAlongMeridian(obj, eccentricities, rightEyeVisualFieldMeridianName, eccUnits, densityUnits);
        
        % Return ISETBio retinal angle for Watson's meridians (specified in visual field of the right eye)
        [isetbioAngle, whichEye, rightEyeRetinalMeridianName] = isetbioRetinalAngleForWatsonMeridian(obj, rightEyeVisualFieldMeridianName);
        
        % Compute 2D cone RF density map for eccentricities specified in the right eye visual space
        [coneRFDensity, spatialSupport, xLabelString, yLabelString, ...
            meridianDensities, densityLabelString, eccUnits, densityUnits] = ...
            compute2DConeRFDensity(obj, eccDegsInREVisualSpace,  theReturnedView, varargin);

        % Compute 2D mRGC RF density map for eccentricities specified in the right eye visual space
        [mRGCRFDensity, spatialSupport, xLabelString, yLabelString, ...
            meridianDensities, densityLabelString, eccUnits, densityUnits] = ...
            compute2DmRGCRFDensity(obj, eccDegsInREVisualSpace,  theReturnedView, varargin);
        
        % Compute 2D cone to mRGC RF ratio map for eccentricities specified in the right eye visual space
        [conesToMRGCratio, spatialSupport,xLabelString, yLabelString, ratioLabel, ...
            meridianConeToMRGratios, eccUnits] = ...
            compute2DConeToMRGCRFRatio(obj, eccDegsInREVisualSpace, eccUnits);
        
        
    end % Public methods
    
    methods (Access=private)

        % Angular interpolation from meridian values to full 360
        val = interpolatedValuesFromMeridianValues(obj, meridianValues, requestedAngles);
        
        % Assert whether the passed eccUnits have valid value
        validateEccUnits(obj,eccUnits);
        
        % Validate densityUnits
        validateDensityUnits(obj,densityUnits);
    
        % Assert whether the passed meridianName has valid value
        validateMeridianName(obj,meridianName);
        
        % Assert whether the passed view has valid value
        validateView(obj, viewName);
    end
    
    % Unit tests
    methods (Static)
        plotLabOBJ = setUpPlotLab();
        plotlabOBJ = unitTestFigure1(varargin);
        plotlabOBJ = unitTestFigure5(varargin);
        plotlabOBJ = unitTestFigure9(varargin);
        plotlabOBJ = unitTestFigure10(varargin);
        plotlabOBJ = unitTestFigure11(varargin);
        plotlabOBJ = unitTestFigure14(varargin);
        plotlabOBJ = unitTestFigureA1(varargin);
        plotlabOBJ = unitTestRetinalSizeToVisualSize(varargin);
        unitTestRFConeDensity2D();
        unitTestRFDensity2D();
        unitTestSmoothGrid();
        samplingVector = generateSamplingVectorFromScatteredXYPositions(scatteredXYPositionsMicrons, nSamples);

    end
    
end % Classdef
    