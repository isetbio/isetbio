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
        % Meridian parameters
        meridianParamsTable = {
            'temporal meridian'  struct('a_k', 0.9851, 'r_2k', 1.058,  'r_ek', 22.14); ...
            'superior meridian'  struct('a_k', 0.9935, 'r_2k', 1.035,  'r_ek', 16.35); ...
            'nasal meridian'     struct('a_k', 0.9729, 'r_2k', 1.084,  'r_ek',  7.633); ...
            'inferior meridian'  struct('a_k', 0.996,  'r_2k', 0.9932, 'r_ek', 12.13);
        }
     
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
       
       % Dictionary with meridian params indexed by meridian name
       meridianParams;
       
       % Struct with default preferences for all figures
       defaultFigurePrefs = struct(...
            'lineWidth', 1.5, ...
            'fontSize', 14, ...
            'fontAngle', 'italic', ...
            'grid', 'on', ...
            'backgroundColor', [1 1 1]);
    end
       
    % Public properties
    properties
        % Struct with default preferences for all figures
        figurePrefs
        
        % the eccentricity support for all figures
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
        
            % Generate figures
            if (generateAllFigures)
                obj.generateAndDockAllFigures();
            end
        end
        
        % --------------------- COMPUTE METHODS ---------------------------
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
        
        % Return cone RF spacing and density at the requested meridian and eccentricities
        [coneRFSpacing, coneRFDensity] = coneRFSpacingAndDensity(obj, eccDegs, meridian, units);
        % --------------------- COMPUTE METHODS ---------------------------
        
        
        % ------------------ FIGURE GENERATION METHODS --------------------
        generateAndDockAllFigures(obj);
        
        % Cone density as a function of eccentricity for all quadrants
        generateFigure1(obj, hFig);
        
        % RF density of all RGCs as a function of eccentricity for all quadrants
        generateFigure5(obj, hFig);
        
        % Fraction of midget to total RGCs RFs as a function of eccentrity
        generateFigure8(obj, hFig);
        
        % RF density of midget RGCs as a function of eccentricity for all quadrants
        generateFigure9(obj, hFig);
        
        % RF spacing of midget RGCs as a function of eccentricity for all quadrants
        generateFigure10(obj, hFig);
        
        % RF spacing of midget RGCs as a function of eccentricity for all quadrants
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
        % ------------------ Figure generation methods --------------------
    end
end

