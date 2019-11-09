function demoWatsonRGCModel
% Demo different ways of using the WatsonRGCModel class
%
% Syntax:
%   demoWatsonRGCModel();
%
% Description:
%  Demo different ways of using the WatsonRGCModel class
%
% References:
%    Watson (2014). 'A formula for human RGC receptive field density as
%    a function of visual field location', JOV (2014), 14(7), 1-17.
%
% History:
%    11/8/19  NPC, ISETBIO Team     Wrote it.

    % Instantiate a WatsonRGCModel object and generate several figures of the Watson 2014 paper
    WatsonRGCCalc = WatsonRGCModel('generateAllFigures', true);
    
    % Compute peak cone density in cones/mm2 and in cones/deg2
    fprintf('peak cone density: %3.2fk cones/mm2\n', WatsonRGCCalc.peakConeDensity('cones per mm2')/1000);
    fprintf('peak cone density: %3.2fk cones/deg2\n', WatsonRGCCalc.peakConeDensity('cones per deg2')/1000);
    
    % Compute peak midget and total RGC RFs
    [peakMidgetRGCRFDensityPerDeg2, peakTotalRGCRFDensityPerDeg2] = WatsonRGCCalc.peakRGCRFDensity('RFs per deg2');
    fprintf('peak midget RGC RF density: %3.2fk RFs/deg2\n', peakMidgetRGCRFDensityPerDeg2/1000);
    fprintf('peak total RGC RF density: %3.2fk RFs/deg2\n', peakTotalRGCRFDensityPerDeg2/1000);
   
    % Compute total RGC RF density along the superior meridian for a number of eccentricities
    eccDegs = [0 0.5 1 1.5];
    meridian = 'superior meridian';
    totalRGCRFDensityPerDeg2 = WatsonRGCCalc.totalRGCRFDensity(eccDegs, meridian, 'RFs per deg2');
    for k = 1:numel(eccDegs)
        fprintf('total RGC RF density at %2.2f degs along %s: %3.2fk RFs/deg2\n', eccDegs(k), meridian, totalRGCRFDensityPerDeg2(k)/1000);
    end


end

