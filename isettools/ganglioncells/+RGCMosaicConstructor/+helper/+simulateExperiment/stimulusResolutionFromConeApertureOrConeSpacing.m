function retinalImageResolutionDegs = stimulusResolutionFromConeApertureOrConeSpacing(...
    theRGCMosaic, targetRGCindices, theFraction, coneApertureOrConeSpacing)

    switch (coneApertureOrConeSpacing)
        case 'cone aperture'
            theData = theRGCMosaic.inputConeMosaic.coneApertureDiametersDegs;
        case 'cone spacing'
            theData = theRGCMosaic.inputConeMosaic.coneRFspacingsDegs;
        otherwise
            error('Invalid coneApertureOrConeSpacing string. Must be either ''cone aperture'' or ''cone spacing''.');
    end

    if (isempty(targetRGCindices))
        retinalImageResolutionDegs = theFraction * min(theData);
    else
        coneIndices = find(theRGCMosaic.rgcRFcenterConeConnectivityMatrix(:,targetRGCindices(1)) > 0.001);
        coneIndices = cat(1, coneIndices, ...
                          find(theRGCMosaic.rgcRFcenterConeConnectivityMatrix(:,targetRGCindices(2)) > 0.001));
        retinalImageResolutionDegs = theFraction * mean(theData(coneIndices)); 
    end

end