function val = peakConeDensity(obj, units)
% Return peak cone density (#cones per either deg^2 or mm^2)
%
% Syntax:
%   WatsonRGCCalc = WatsonRGCModel();
%   peakConeDensityPerMM2 = WatsonRGCCalc.peakConeDensity('cones per mm2')
%   peakConeDensityPerDeg2 = WatsonRGCCalc.peakConeDensity('cones per deg2')
%
% Description:
%   Method to return the peak cone density (#of cones/area, with area 
%   specified either in deg^2 or mm^2.
%
% Inputs:
%    obj                       - The WatsonRGCModel object
%    units                     - Retinal area units, either 'cones per mm2'
%                                or 'cones per deg2'
% Outputs:
%    val                       - Peak cone density in specified units
%
% References:
%    Watson (2014). 'A formula for human RGC receptive field density as
%    a function of visual field location', JOV (2014), 14(7), 1-17.
%
% History:
%    11/8/19  NPC, ISETBIO Team     Wrote it.

    % Retrieve peak cone density in cones/deg2
    coneDensityPerDeg2 = obj.dc0;
    
    switch (units)
        case 'Cones per mm2'
            % Compute mmSquaredPerDegSquared conversion factor at zero eccentricity
            eccDegs = 0;
            mmSquaredPerDegSquared = obj.alpha(eccDegs);
            % Compute degSquaredPerMMSquared conversion factor
            deg2PerMM2 = 1/mmSquaredPerDegSquared;
            % Compute cones per mm2
            val = coneDensityPerDeg2 * deg2PerMM2;
        case 'Cones per deg2'
            val = coneDensityPerDeg2;
        otherwise
            error('units must be either ''Cones per mm2'' or ''Cones per deg2''. Instead it was ''%s''.', units);
    end
end