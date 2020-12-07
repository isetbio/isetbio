classdef constants
% Constants defined in Watson (2014): ''A formula for human RGC receptive field density as a function of visual field location'

    properties (Constant)

        % Peak cone density (cones/deg^2) at 0 deg eccentricity (page 3, dc(0), Also in Appendix 4)
        dc0 = 14804.6;
        
        % Fraction of all ganglion cells that are midget at 0 deg eccentricity
        f0 = 1/1.12;
        
        % Fraction of all ganglion cells that are midget as a function of
        % eccentricity in degrees (Equation 7)
        midgetRGCFractionEccVariation = @(eccDegs) RGCmodels.Watson.constants.f0 ./ (1 + eccDegs/41.03);

        % Ratio of cone aperture to cone diameter (set to 0.7 in isetbio)
        coneApertureToDiameterRatio = 0.7;
        
        % Label for left eye
        leftEye = 'left eye';
        
        % Label for right eye
        rightEye = 'right eye';
        
        % Label for nasal meridian
        nasalMeridian    = 'nasal meridian';
        
        % Label for temporal meridian
        temporalMeridian = 'temporal meridian';
        
        % Label for inferior meridian
        inferiorMeridian = 'inferior meridian';
        
        % Label for superior meridian
        superiorMeridian = 'superior meridian';
        
        % Right eye visual field meridian enumeration (mathing Watson's paper)
        indexedMeridians = {...
            RGCmodels.Watson.constants.temporalMeridian ...
            RGCmodels.Watson.constants.superiorMeridian ...
            RGCmodels.Watson.constants.nasalMeridian ...
            RGCmodels.Watson.constants.inferiorMeridian ...
            };
        
        % Right eye visual field meridian angles (matching Watson's paper)
        indexedMeridianAngles = [0 90 180 270];
        
        % Meridian params for function modeling the variation of RGC density 
        % with eccentrcitity (Table 1)
        meridianParamsTable = containers.Map(...
            RGCmodels.Watson.constants.indexedMeridians, ...
            {struct('a_k', 0.9851, 'r_2k', 1.058,  'r_ek', 22.14), ... 
             struct('a_k', 0.9935, 'r_2k', 1.035,  'r_ek', 16.35), ...
             struct('a_k', 0.9729, 'r_2k', 1.084,  'r_ek',  7.633), ... 
             struct('a_k', 0.996,  'r_2k', 0.9932, 'r_ek', 12.13)});
         
        % Function modeling the variation of RGC density with eccentricity
        % (Equation 4)
        totalRGCRFDensityEccVariation = @(eccDegs, params) ...
              params(1) * (1 + eccDegs/params(2)).^(-2) ...
            + (1-params(1)) * exp(-eccDegs/params(3));
        
        % Right eye visual field meridian colors (matching Watson's paper figures)
        meridianColors = [...
            [1 0 0]; ...
            [0 0 1]; ...
            [0 0.7 0]; ...
            [0.3 0.3 0.3] ...
            ];
    end
    
end
