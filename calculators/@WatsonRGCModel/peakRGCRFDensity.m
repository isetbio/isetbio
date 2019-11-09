function [peakMidgetRGCRFDensity, peakTotalRGCRFDensity] = peakRGCRFDensity(obj, units)
% Return the peak midget RGC and peak total RGC RF densities (#RFs per either deg^2 or mm^2)
%
% Syntax:
%   WatsonRGCCalc = WatsonRGCModel();
%   [peakMidgetRGCRFDensityPerMM2, peakTotalRGCRFDensityPerMM2] = WatsonRGCCalc.peakRGCRFDensity('RFs per mm2')
%   [peakMidgetRGCRFDensityPerDeg2, peakTotalRGCRFDensityPerDeg2] = WatsonRGCCalc.peakRGCRFDensity('RFs per deg2')
%
% Description:
%   Method to return the peak midget RGC and peak total RGC RF densities
%   (#of RFs/area, with area specified either in deg^2 or mm^2.
%
% Inputs:
%    obj                       - The WatsonRGCModel object
%    units                     - Retinal area units, either 'RFs per mm2'
%                                or 'RFs per deg2'
% Outputs:
%    peakMidgetRGCRFDensity    - Peak midget RGC RF density in specified units
%    peakTotalRGCRFDensity     - Peak total RGC RF density in specified units
% 
% References:
%    Watson (2014). 'A formula for human RGC receptive field density as
%    a function of visual field location', JOV (2014), 14(7), 1-17.
%
% History:
%    11/8/19  NPC, ISETBIO Team     Wrote it.

    % Retrieve cone density
    if (contains(units, 'per mm2'))
        peakConeDensity = obj.peakConeDensity('cones per mm2');
    elseif (contains(units, 'per deg2'))
        peakConeDensity = obj.peakConeDensity('cones per deg2');
    else
        error('units must be either ''RFs per mm2'' or ''RFs per deg2''.');
    end
    
    % The foveal density of midget RGC RFs is twice the cone density
    % because in the fovea each cone connects to exactly 2 midget RGCs (one
    % ON and one OFF). This is Equation (1) in the Watson (2014) paper.
    peakMidgetRGCRFDensity = 2 * peakConeDensity;
    
    % Retrieve percentage of total ganglion cells that are midget at 0 eccentricity
    percentageOfRGCsThatAreMidgetsAtZeroEccentricity = obj.f0;
    % This is Equation (2) in the Watson (2014) paper.
    peakTotalRGCRFDensity = 1/percentageOfRGCsThatAreMidgetsAtZeroEccentricity * peakMidgetRGCRFDensity;
end