function val = midgetRGCRFDensity(obj, eccDegs, meridian, units)
% Return midget RGC receptive field density at the requested meridian and eccentricities
%
% Syntax:
%   WatsonRGCCalc = WatsonRGCModel();
%   eccDegs = 0:0.1:10;
%   meridian = 'superior meridian';
%   midgetRGCRFDensityPerMM2 = WatsonRGCCalc.midgetRGCRFDensity(eccDegs, meridian, 'RFs per mm2')
%   midgetRGCRFDensityPerDeg2 = WatsonRGCCalc.midgetRGCRFDensity(eccDegs, meridian, 'RFs per deg2')
%
% Description:
%   Method to return the midget RGC receptive field density (ON and OFF submosaics combined)
%   as a function of the requested meridian and eccentricities (#of RFs/area, with area specified either in deg^2 or mm^2.
%
% Inputs:
%    obj                       - The WatsonRGCModel object
%    eccDegs                   - Eccentricities at which to compute RF densities
%    meridian                  - Meridian for which to compute RF densities
%    units                     - Retinal area units, either 'RFs per mm2'
%                                or 'RFs per deg2'
% Outputs:
%    val                       - Midget RGC receptive field densities at the requested eccentricities
% 
% References:
%    Watson (2014). 'A formula for human RGC receptive field density as
%    a function of visual field location', JOV (2014), 14(7), 1-17.
%
% History:
%    11/11/19  NPC, ISETBIO Team     Wrote it.
    
    % Compute total RGC receptive field  density at the requested units
    totalRGCRFdensity = obj.totalRGCRFDensity(eccDegs, meridian, units);
    midgetRGCFraction = obj.midgetRGCFraction(eccDegs);
    
    % This is equation (8) in the Watson (2014) paper.
    val = midgetRGCFraction .* totalRGCRFdensity;
end