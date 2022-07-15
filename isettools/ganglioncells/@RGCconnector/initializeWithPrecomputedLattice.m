function RGCRFposMicrons = initializeWithPrecomputedLattice(obj, modelRGC)

    assert(ismember(modelRGC, obj.validRGClatticeModels), 'Unknown modelRGC: ''%s''.', modelRGC);

    % Generate RGCRF positions that lie within the limits of the input cone mosaic
    allConePositions = obj.inputConeMosaic.coneRFpositionsMicrons;
    minConePosX = min(allConePositions(:,1));
    minConePosY = min(allConePositions(:,2));
    maxConePosX = max(allConePositions(:,1));
    maxConePosY = max(allConePositions(:,2));

    eccMicrons = [mean(allConePositions(:,1)) mean(allConePositions(:,2))];
    sizeMicrons = max([maxConePosX-minConePosX maxConePosY-minConePosY]);

    % Import positions from a pre-generated RGC RF lattice
    
    switch (modelRGC)
        case 'Watson-midgetRGC'
            sourceLatticeSizeDegs = 58;
            customDegsToMMsConversionFunction = @(x) RGCmodels.Watson.convert.rhoDegsToMMs(x);
    
            RGCRFposMicrons = retinalattice.import.finalMRGCPositions(...
                sourceLatticeSizeDegs, ...
                eccMicrons, ...
                sizeMicrons, ...
                obj.inputConeMosaic.whichEye, ...
                customDegsToMMsConversionFunction); 
        otherwise
            error('No data for importing ''%s''-class lattice', modelRGC);
    end
end