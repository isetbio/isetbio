function [coneRFSpacing, coneRFDensity] = coneRFSpacingAndDensity(obj, eccDegs, meridian, units)
% Return cone receptive field spacing and density at the requested meridian and eccentricities
%
% Syntax:
%   WatsonRGCCalc = WatsonRGCModel();
%   eccDegs = 0:0.1:10;
%   meridian = 'superior meridian';
%   [coneSpacingMM, coneDensityPerMM2] = WatsonRGCCalc.coneRFSpacingAndDensity(eccDegs, meridian, 'mm')
%   [coneSpacingDeg, coneDensityPerDeg2] = WatsonRGCCalc.coneRFSpacingAndDensity(eccDegs, meridian, 'deg')
%
% Description:
%   Method to return cone spacing as a function of the requested meridian and 
%   eccentricities (with spacing specified either in deg or mm).
%
% Inputs:
%    obj                       - The WatsonRGCModel object
%    eccDegs                   - Eccentricities at which to compute RF densities
%    meridian                  - Meridian for which to compute RF densities
%    units                     - Retinal area units, either 'RFs per mm2'
%                                or 'RFs per deg2'
% Outputs:
%    val                       - Cone spacing at the requested eccentricities
% 
% References:
%    Watson (2014). 'A formula for human RGC receptive field density as
%    a function of visual field location', JOV (2014), 14(7), 1-17.
%
% History:
%    11/11/19  NPC, ISETBIO Team     Wrote it.
    
    % Load the Curcio '1990 cone spacing data
    switch (meridian)
        case 'temporal meridian'
            angle = 180;
        case 'superior meridian'
            angle = 90;
        case 'nasal meridian'
            angle = 0;
        case 'inferior meridian'
            angle = 270;
        otherwise
            fprintf('Valid meridian names are: %s\n', keys(obj.meridianParams));
            error('Invalid meridian name: ''%s''.', meridian);
    end
    
    % Convert ecc in degs to ecc in MMs
    eccMM = obj.rhoDegsToMMs(eccDegs);
    
    % Call the isetbio function coneSizeReadData to read-in the Curcio '1990
    % cone spacing/density data
    [spacingMeters, apertureMeters, densityConesPerMM2] = coneSizeReadData('eccentricity', eccMM, ...
                                        'angle', angle*ones(1,numel(eccMM)), ...
                                        'eccentricityUnits', 'mm', ...
                                        'useParfor', false);
    
    switch (units)
        case 'deg'
            spacingMM = spacingMeters * 1e-3;
            % Convert cone spacing in mm to cone spacing in degs at all eccentricities
            coneRFSpacing = obj.rhoMMsToDegs(spacingMM+eccMM)-obj.rhoMMsToDegs(eccMM); 
            
            % Convert cone density from per mm2 to per deg2
            % Compute mmSquaredPerDegSquared conversion factor for the eccentricities in mm
            mmSquaredPerDegSquared = obj.alpha(eccMM);
            coneRFDensity = densityConesPerMM2 .* mmSquaredPerDegSquared;
        case 'mm'
            coneRFSpacing = spacingMeters * 1e-3;
            coneRFDensity = densityConesPerMM2;
        otherwise
            error('Units must be either ''deg'' or ''mm''.');
    end
    
end