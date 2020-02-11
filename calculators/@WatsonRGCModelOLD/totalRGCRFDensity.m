function val = totalRGCRFDensity(obj, eccDegs, meridian, units)
% Return total RGC receptive field  density at the requested meridian and eccentricities
%
% Syntax:
%   WatsonRGCCalc = WatsonRGCModel();
%   eccDegs = 0:0.1:10;
%   meridian = 'superior meridian';
%   totalRGCRFDensityPerMM2 = WatsonRGCCalc.totalRGCRFDensity(eccDegs, meridian, 'RFs per mm2')
%   totalRGCRFDensityPerDeg2 = WatsonRGCCalc.totalRGCRFDensity(eccDegs, meridian, 'RFs per deg2')
%
% Description:
%   Method to return the total RGC receptive field density as a function of the requested
%   meridian and eccentricities (#of RFs/area, with area specified either in deg^2 or mm^2.
%
% Inputs:
%    obj                       - The WatsonRGCModel object
%    eccDegs                   - Eccentricities at which to compute RF densities
%    meridian                  - Meridian for which to compute RF densities
%    units                     - Retinal area units, either 'RFs per mm2'
%                                or 'RFs per deg2'
% Outputs:
%    val                       - Total RGC receptive field densities at the requested eccentricities
% 
% References:
%    Watson (2014). 'A formula for human RGC receptive field density as
%    a function of visual field location', JOV (2014), 14(7), 1-17.
%
% History:
%    11/8/19  NPC, ISETBIO Team     Wrote it.

    % Retrieve peak total RGC receptive field density at the requested units
    [~, peakTotalRGCRFDensity] = obj.peakRGCRFDensity(units);
    
    meridianParams = obj.meridianParams(meridian);
    
    % This is equation (4) in the Watson (2014) paper.
    eccVariation = ...
        meridianParams.a_k * (1+eccDegs/meridianParams.r_2k).^(-2) + ...
        (1-meridianParams.a_k) * exp(-eccDegs/meridianParams.r_ek);

    val = peakTotalRGCRFDensity * eccVariation;
end