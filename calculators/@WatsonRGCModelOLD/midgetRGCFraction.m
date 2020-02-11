function val = midgetRGCFraction(obj, eccDegs)
% Return fraction of total RGCs that are midgets
%
% Syntax:
%   WatsonRGCCalc = WatsonRGCModel();
%   eccDegs = 0:0.1:10;
%   midgetRGCf = WatsonRGCCalc.midgetRGCFraction(eccDegs)
%
% Description:
%   Method to return the fraction of total RGCs that are midgets at a range
%   of eccentricities.
%
% Inputs:
%    obj                       - The WatsonRGCModel object
%    eccDegs                   - Eccentricities at which to compute RF densities
%
% Outputs:
%    val                       - Fraction of total RGCs that are midgets at the requested eccentricities
% 
% References:
%    Watson (2014). 'A formula for human RGC receptive field density as
%    a function of visual field location', JOV (2014), 14(7), 1-17.
%
% History:
%    11/11/19  NPC, ISETBIO Team     Wrote it.

    % Retrieve percentage of total ganglion cells that are midget at 0 eccentricity
    percentageOfMidgetRGCsAtZeroEccentricity = obj.f0;
    eccMDegs = 41.03;
    
    % This is equation (7) in the Watson (2014) paper.
    val = percentageOfMidgetRGCsAtZeroEccentricity ./ (1 + eccDegs/eccMDegs);
end