function [mRGCDensityPerMM2, maxDensity] = mRGCRFDensityFunctionFull(rfPositions, whichEye)

    % Instantiate a WatsonRGCModel
    WatsonOBJ = WatsonRGCModel();
    
    rfPositionsMM = 1e-3*rfPositions;
    eccUnits = 'mm';
    densityUnits = 'mm^2';
    
    [~,mRGCDensityPerMM2] = WatsonOBJ.mRGCRFSpacingAndDensityAtRetinalPositions(rfPositionsMM, whichEye, ...
        eccUnits, densityUnits, 'adjustForISETBioConeDensity', true, ...
        'subtype', 'ON'); 
    mRGCDensityPerMM2 = mRGCDensityPerMM2';
    
    [~,maxDensity] = WatsonOBJ.mRGCRFSpacingAndDensityAtRetinalPositions([0 0], whichEye, ...
         eccUnits, densityUnits, 'adjustForISETBioConeDensity', true, ...
         'subtype', 'ON');
end