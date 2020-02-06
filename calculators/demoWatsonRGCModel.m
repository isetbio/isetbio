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

    % Instantiate a WatsonRGCModel object. Set the 'generateAllFigures' flag 
    % to true to generate several figures of the Watson 2014 paper
    WatsonRGCCalc = WatsonRGCModel('generateAllFigures', false);
    
    % How to sample visual space
    eccMinArcMin = 0.2;
    eccMaxDegs = 90;
    eccSamples = 200;
    WatsonRGCCalc.eccDegs = logspace(log10(eccMinArcMin/60), log10(eccMaxDegs), eccSamples);
    
    % Generate the meridian conventions figure
    WatsonRGCCalc.generateMeridianConventionsFigure();
    
    % Plot cone density as a function of ecc for the left eye (by default, meridians are labeled in the right eye visual space)
    WatsonRGCCalc.generateFigure1(figure(10), 'whichEye', 'left', ...
        'eccentricityInMMInsteadOfDegs', false);
    
    % Plot cone density as a function of ecc for the left eye with meridians labeled in the eye's own retinal space
    WatsonRGCCalc.generateFigure1(figure(11), 'whichEye', 'left', ...
        'eccentricityInMMInsteadOfDegs', false, 'retinalMeridiansLegendsInsteadOfVisualSpaceMeridians', true);
    
    % Plot cone density as a function of ecc for the left eye with meridians labeled in the eye's own retinal space 
    % and in retinal mm instead of visual degs
    WatsonRGCCalc.generateFigure1(figure(12), 'whichEye', 'left', ...
        'eccentricityInMMInsteadOfDegs', true, 'retinalMeridiansLegendsInsteadOfVisualSpaceMeridians', true);
    
    % Plot cone density as a function of ecc for the right eye (by default, meridians are labeled in the right eye visual space)
    WatsonRGCCalc.generateFigure1(figure(20), 'whichEye', 'right', ...
        'eccentricityInMMInsteadOfDegs', false);
    
    % Plot cone density as a function of ecc for the lright eye with meridians labeled in the eye's own retinal space
    WatsonRGCCalc.generateFigure1(figure(21), 'whichEye', 'right', ...
        'eccentricityInMMInsteadOfDegs', false, 'retinalMeridiansLegendsInsteadOfVisualSpaceMeridians', true);
    
    % Plot cone density as a function of ecc for the left eye with meridians labeled in the eye's own retinal space 
    % and in retinal mm instead of visual degs
    WatsonRGCCalc.generateFigure1(figure(22), 'whichEye', 'right', ...
        'eccentricityInMMInsteadOfDegs', true, 'retinalMeridiansLegendsInsteadOfVisualSpaceMeridians', true);
    
end

function oldTests()
    % Compute peak cone density in cones/mm2 and in cones/deg2
    fprintf('peak cone density: %3.2fk cones/mm2\n', WatsonRGCCalc.peakConeDensity('Cones per mm2')/1000);
    fprintf('peak cone density: %3.2fk cones/deg2\n', WatsonRGCCalc.peakConeDensity('Cones per deg2')/1000);
    
    % Compute peak midget and total RGC RFs
    [peakMidgetRGCRFDensityPerDeg2, peakTotalRGCRFDensityPerDeg2] = WatsonRGCCalc.peakRGCRFDensity('RFs per deg2');
    fprintf('peak midget RGC RF density: %3.2fk RFs/deg2\n', peakMidgetRGCRFDensityPerDeg2/1000);
    fprintf('peak total RGC RF density: %3.2fk RFs/deg2\n', peakTotalRGCRFDensityPerDeg2/1000);
   
    % Compute total RGC RF density along the superior meridian for a number of eccentricities
    eccDegs = [0.2 0.5 1 1.5];
    meridian = 'superior meridian';
    totalRGCRFDensityPerDeg2 = WatsonRGCCalc.totalRGCRFDensity(eccDegs, meridian, 'RFs per deg2');
    for k = 1:numel(eccDegs)
        fprintf('total RGC RF density at %2.2f degs along %s: %3.2fk RFs/deg2\n', eccDegs(k), meridian, totalRGCRFDensityPerDeg2(k)/1000);
    end

    % Compute ratio of midget RGC RFs to cones for a number of eccentricities
    eccDegs = [0.2 0.5 1 2 3 4 10];
    midgetRGCRFDensity = WatsonRGCCalc.midgetRGCRFDensity(eccDegs, meridian, 'RFs per deg2');
    [coneSpacingDegs, coneDensity] = WatsonRGCCalc.coneRFSpacingAndDensity(eccDegs, meridian, 'Cones per deg2');
    midgetRGCRFtoConeRatio = midgetRGCRFDensity ./ coneDensity;
    for k = 1:numel(eccDegs)
        fprintf('midget RGC RF (ON+OFF) -to- cone ratio at %2.1f degs along %s: %3.2f (cone spacing: %2.2f arc min)\n', eccDegs(k), meridian, midgetRGCRFtoConeRatio(k), coneSpacingDegs(k)*60);
    end
    
    % Compute the average number of cones per midget RGC RF ratio for a number of eccentricities
    % The density of the ON or OFF midget RGC mosaic is half the density of
    % the total midgetRGC mosaic assuming equal numerosities of the ON and
    % the OFF submosaics
    midgetRGCRFSinglePolarityDensity = midgetRGCRFDensity/2;
    % Compute the average # of cones per midget RGC by dividing the cone
    % density by the midgetRGC RF density for the ON/OFF submosaic
    averageNumberOfConesPerMidgetRGCcenter = coneDensity ./ midgetRGCRFSinglePolarityDensity;
    for k = 1:numel(eccDegs)
        fprintf('average number of cones per midget RGC RF (ON or OFF) center at %2.1f degs along %s: %2.2f\n', eccDegs(k), meridian, averageNumberOfConesPerMidgetRGCcenter(k));
    end
    
end

