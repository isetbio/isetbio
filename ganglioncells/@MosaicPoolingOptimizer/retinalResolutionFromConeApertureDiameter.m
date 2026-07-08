function retinalImageResolutionDegs = retinalResolutionFromConeApertureDiameter(theRGCMosaic, targetRGCindices)

    if (isempty(targetRGCindices))
        retinalImageResolutionDegs = min(theRGCMosaic.inputConeMosaic.coneApertureDiametersDegs)/9;
    else
        % Find cone aperture size of the input cones
        if (~isempty(theRGCMosaic.rgcRFcenterConeConnectivityMatrix))
            coneIndices = find(theRGCMosaic.rgcRFcenterConeConnectivityMatrix(:,targetRGCindices(1)) > 0.001);
            coneIndices = cat(1, coneIndices, ...
                          find(theRGCMosaic.rgcRFcenterConeConnectivityMatrix(:,targetRGCindices(2)) > 0.001));
        else
            coneIndices = find(theRGCMosaic.rgcRFcenterConePoolingMatrix(:,targetRGCindices(1)) > 0.001);
            coneIndices = cat(1, coneIndices, ...
                          find(theRGCMosaic.rgcRFcenterConePoolingMatrix(:,targetRGCindices(2)) > 0.001));
        end
        
        retinalImageResolutionDegs = mean(theRGCMosaic.inputConeMosaic.coneApertureDiametersDegs(coneIndices))/13; 
    end

end