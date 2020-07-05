function runPhaseX(runParams)

    % Intermediate files directory
    saveDir = strrep(fileparts(which(mfilename())), 'processing', 'responseFiles');
    
    % Compute cone mosaic responses
    recomputeConeMosaicResponses = ~true;
    recomputeNullResponses = ~true;
    recomputeMosaicsAndOptics = ~true;
    
    % Load/Recompute connected mosaics and the optics
    [theConeMosaic, theMidgetRGCmosaic, theOptics] = mosaicsAndOpticsForEccentricity(runParams, recomputeMosaicsAndOptics, saveDir);

    % Stimulation parameters
    LMScontrast = [0.1 0.0 0.0];
    minSF = 0.1;
    maxSF = 60;
    spatialFrequenciesCPD = logspace(log10(minSF), log10(maxSF),12);
    
    stimulusFOVdegs = 2.0;
    minPixelsPerCycle = 10;
    stimulusPixelsNum = maxSF*stimulusFOVdegs*minPixelsPerCycle;
    temporalFrequency = 2.0;
    stimDurationSeconds = 1.0;
    instancesNum = 2;
    
    
    
    % Signal to the RGCs
    rgcInputSignal = 'isomerizations';
    %rgcInputSignal = 'photocurrents';
    
    if (recomputeConeMosaicResponses)
        computeConeResponses(runParams, theConeMosaic, theOptics, ...
            recomputeNullResponses, ...
            instancesNum, stimDurationSeconds, stimulusFOVdegs, stimulusPixelsNum, ...
            spatialFrequenciesCPD, temporalFrequency, LMScontrast, ...
            saveDir);
    else
        computeRGCresponses(runParams, theConeMosaic, theMidgetRGCmosaic, ...
            rgcInputSignal, spatialFrequenciesCPD, LMScontrast, ...
            saveDir);
    end
end