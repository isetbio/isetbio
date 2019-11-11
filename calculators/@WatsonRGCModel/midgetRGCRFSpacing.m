function val = midgetRGCRFSpacing(obj, eccDegs, meridian, type, units)
% Return midget RGC receptive field spacing at the requested meridian and eccentricities for a given type (ON/OFF or both ON+OFF).
%
% Syntax:
%   WatsonRGCCalc = WatsonRGCModel();
%   eccDegs = 0:0.1:10;
%   meridian = 'superior meridian';
%   type = 'ON';
%   type = 'single polarity';
%   midgetRGCRFSpacingMM = WatsonRGCCalc.midgetRGCRFSpacing(eccDegs, meridian, type, 'mm')
%   midgetRGCRFSpacingDeg = WatsonRGCCalc.midgetRGCRFSpacing(eccDegs, meridian, type, 'deg')
%
% Description:
%   Method to return the midget RGC receptive field spacing as a function of the requested
%   meridian and eccentricities for a given type ((ON/OFF or both ON+OFF) 
%   (with spacing specified either in deg or mm).
%
% Inputs:
%    obj                       - The WatsonRGCModel object
%    eccDegs                   - Eccentricities at which to compute RF densities
%    meridian                  - Meridian for which to compute RF densities
%    units                     - Retinal spacing units, either 'mm' or 'deg''
%
% Outputs:
%    val                       - Midget RGC receptive field spacing at the requested eccentricities
% 
% References:
%    Watson (2014). 'A formula for human RGC receptive field density as
%    a function of visual field location', JOV (2014), 14(7), 1-17.
%
% History:
%    11/11/19  NPC, ISETBIO Team     Wrote it.

    % Compute total midget RGC receptive field density at the requested units
    switch (units)
        case 'mm'
            densityUnits = 'RFs per mm2';
        case 'deg'
            densityUnits = 'RFs per deg2';
        otherwise
            error('Units for spacing must be either ''deg'' or ''mm''.');
    end
    
    if (~ismember(type, {'single polarity', 'both polarities'}))
        error('RGC polarity for spacing must be either ''single polarity'' (ON/OFF) or ''both polarities'' (ON+OFF).');
    end
    
    % Compute midget RGC receptive field densities.
    midgetRGCRFDensity = obj.midgetRGCRFDensity(eccDegs, meridian, densityUnits);
    
    % Half the RF density if density if requested for a single polarity
    % channel (ON or OFF). We are ignoring assymetries in ON/OFF RGC numbers.
    if (strcmp(type, 'single polarity'))
        midgetRGCRFDensity = midgetRGCRFDensity/2;
    end
    
    % This is equation (9) in the Watson (2014) paper.
    val = sqrt(2.0./(sqrt(3.0)*midgetRGCRFDensity));
end