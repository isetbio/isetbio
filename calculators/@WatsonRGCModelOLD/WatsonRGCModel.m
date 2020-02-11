classdef WatsonRGCModel
    % Create a WatsonRGCModel
    %
    %
    % Syntax:
    %   cMosaic = WatsonRGCModel('generateAllFigures', false);
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
        meridianParamsTable = {
            'temporal meridian'  struct('a_k', 0.9851, 'r_2k', 1.058,  'r_ek', 22.14); ... 
            'superior meridian'  struct('a_k', 0.996,  'r_2k', 0.9932, 'r_ek', 12.13); ...
            'nasal meridian'     struct('a_k', 0.9729, 'r_2k', 1.084,  'r_ek',  7.633); ... 
            'inferior meridian'  struct('a_k', 0.9935, 'r_2k', 1.035,  'r_ek', 16.35) ... 
        };
     
        
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
        
        % Conversion factor, rho, of retinal distance deg->mm as as a function of eccentricity in
        % degs (Equation A5)
        rhoDegsToMMs = @(eccDegs) ...
            0.268         * eccDegs + ...
            0.0003427     * eccDegs .^2 + ...
           -8.3309 * 1e-6 * eccDegs .^3;
        
        % Conversion factor, rho, of retinal distance mm->deg as as a function of eccentricity in
        % degs (Equation A6)
        rhoMMsToDegs = @(eccMM) ...
            3.556     * eccMM + ...
            0.05993   * eccMM .^2 + ...
           -0.007358  * eccMM .^3 + ...
            0.0003027 * eccMM .^4;
       
       
        % Conversion factor, alpha, of retinal area mm^2 -> deg^2 as a function of eccentricity in
        % degs (Equation A7)
        alpha = @(eccDegs) 0.0752 + ...
                    5.846 * 1e-5 * eccDegs    + ...
                   -1.064 * 1e-5 * eccDegs.^2 + ...
                    4.116 * 1e-8 * eccDegs.^3;
                
        
    end
   
    % Constant properties related to figure generation
    properties (Constant)
        paperTitleFull = 'Watson (2014): ''A formula for human RGC receptive field density as a function of visual field location'' ';
        paperTitleShort = 'Watson (2014) RGC model';
    end

    % Public properties (read-only)
    properties (SetAccess = private)
       % Dictionary with various acronyms of the the Watson (2014) paper and their meaning
       glossary;
       
       % Enumerated meridian names:  temporal, superior, nasal, inferior
       enumeratedMeridianNames;
        
       % Dictionary with meridian params indexed by meridian name
       meridianParams;
       
       % Dictionary with meridian colors indexed by meridian name
       meridianColors;
       
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
            
            % Create enumerated meridian names
            obj.enumeratedMeridianNames = cell(size(obj.meridianParams,1),1);
            for k = 1:size(obj.meridianParams,1)
                obj.enumeratedMeridianNames{k} = WatsonRGCModel.meridianParamsTable{k,1};
            end
    
            % Create dictionary with meridian colors
            obj.meridianColors = containers.Map();
            obj.meridianColors('temporal meridian') = [1.0 0.0 0.0];
            obj.meridianColors('superior meridian') = [0.0 0.0 1.0];
            obj.meridianColors('nasal meridian')    = [0.0 0.8 0.0];
            obj.meridianColors('inferior meridian') = [0.2 0.2 0.2];
            
            % Generate figures
            if (generateAllFigures)
                obj.generateAndDockAllFigures();
            end
        end
        
        % --------------------- PARAMETER VALIDATION METHODS --------------
        function validateMeridianName(obj, meridianName)
            assert(ismember(meridianName, obj.enumeratedMeridianNames), ...
                sprintf('''%s'' is not a valid meridian name.', meridianName));
        end
        
        function validateMeridianSpace(obj, meridianSpace)
            assert(ismember(meridianSpace, {'retinal', 'visual'}), ...
                sprintf('''%s'' is not a valid meridian space. Choose ''visual'' or ''retinal''. ', meridianSpace));
        end
        
        function validateEye(obj, whichEye)
            assert(ismember(whichEye, {'left', 'right'}), ...
                sprintf('''%s'' is not a valid eye. Choose ''left'' or ''right''. ', whichEye));
        end
        
        % --------------------- COMPUTE METHODS ---------------------------
        
        % Return meridian angles for desired meridian, meridian space and eye
        angles = meridianAngles(obj, meridianName, meridianSpace, whichEye);
        
        %Compute ratios of midget RGC receptive fields to cones at the requested (x,y) retinal eccentricities (degs)
        [ratios, coneRFDensities, mRGCRFDensities] = ratioOfMidgetRGCsToCones(obj, eccXYposDegs, whichEye);
        
        % Convert retinal area from deg^2 to mm^2 for a given eccentricity
        val = mmSquaredToDegSquared(obj, mmSquared, eccDegs);
        
        % Convert retinal area from deg^2 to mm^2 for a given eccentricity
        val = degSquaredToMMSquared(obj, degSquared, eccDegs);
        
        % Return peak cone density (#cones per either deg^2 or mm^2)
        val = peakConeDensity(obj, units);
        
        % Return peak midget and peak total RGC receptive field density 
        [peakMidgetRGCRFDensity, peakRGCRFDensity] = peakRGCRFDensity(obj, units);
        
        % Return total RGC receptive field  density at the requested meridian and eccentricities 
        val = totalRGCRFDensity(obj, eccDegs, meridian, units);
        
        % Return fraction of total RGCs that are midgets
        val = midgetRGCFraction(obj, eccDegs);
        
        % Return midgetRGC receptive field  density at the requested meridian and eccentricities
        val = midgetRGCRFDensity(obj, eccDegs, meridian, units);
        
        % Return midgetRGC receptive field  spacing at the requested meridian and
        % eccentricities for a given RGC type: 'singlePolarity' (ON/OFF) or
        % 'bothPolarities' (ON+OFF).
        val = midgetRGCRFSpacing(obj, eccDegs, meridian, units, type);
        
        % Return cone RF spacing, cone RF density and meridian angle for the 
        % requested meridian and eccentricities
        [coneRFSpacing, coneRFDensity, meridianAngle] = coneRFSpacingAndDensity(obj, ...
            ecc, meridian, meridianSpace, whichEye, eccUnits, returnUnits);
        % --------------------- COMPUTE METHODS ---------------------------
        
        
        % ------------------ FIGURE GENERATION METHODS --------------------
        generateAndDockAllFigures(obj);
        
        % Meridian conventions figure
        generateMeridianConventionsFigure(obj);
        
        % Cone density (cones/deg2) as a function of eccentricity in degs for all quadrants
        generateFigure1(obj, hFig, varargin);
        
        % RF density of all RGCs (RFs/deg2) as a function of eccentricity for all quadrants
        generateFigure5(obj, hFig);
        
        % Fraction of midget to total RGCs RFs as a function of eccentrity
        generateFigure8(obj, hFig);
        
        % RF density of midget RGCs (RFs/deg2) as a function of eccentricity for all quadrants
        generateFigure9(obj, hFig);
        
        % RF spacing of midget RGCs (degs) as a function of eccentricity for all quadrants
        generateFigure10(obj, hFig);
        
        % RF spacing of midget RGCs (degs) as a function of eccentricity for all quadrants
        generateFigure11(obj, hFig);
        
        % Ratio of midget RGCs to cones as a function of eccentricity for all quadrants
        generateFigure14(obj, hFig);
                
        % Relationhip between retinal distance from the optic axis in mm and degs as a
        % function of eccentricity
        generateFigureA1(obj, hFig);
        
        % Ratio of area in mm^2 to deg^2 as a function of eccentricity
        generateFigureA2(obj, hFig);
        
        % Method to generate RGCdensity at four quadrants as a function of eccentricity
        generateRGCRFDensityPlot(obj, RGCRFDensityFunctionHandle, eccDegs);
        
        % Method to generate RGCspacing at four quadrants as a function of
        % eccentricity for a given cell type (ON/OFF or both ON+OFF)
        generateRGCRFSpacingPlot(obj, RGCRFSpacingFunctionHandle, eccDegs, type)
        
        % Correct meridian color & legend based on which eye we are showing and whether we
        % are labeling meridians in retinal or visual space
        [theLegend, theColor] = correctLegendAndColor(obj,theLegend, theColor, meridianName, displayRetinalMeridiansLegends, whichEye)
        
        % ------------------ Figure generation methods --------------------
    end
end

