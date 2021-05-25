function stimSizeDegs = stimulusSizeAtLowestSpatialFrequency(app)

    if (strcmp(app.csfParams.constantParameter, 'constant size'))
        stimSizeDegs = app.stimParams.sizeDegs;
    else
        spatialPeriodDegs = 1.0/app.csfParams.spatialFrequencyMin;
        stimSizeDegs = spatialPeriodDegs * app.csfParams.numberOfConstantCycles*6;
    end
end