function theRF = generateRGCRFcenterSubregionFromConnectivityMatrix(connectivityVectorForRGC, conePositions, coneSpacings,  X,Y)
    
    theRF = [];
    connectedConeIDs = find(connectivityVectorForRGC>0);
    flatTopZ = 0.4;
    
    for k = 1:numel(connectedConeIDs)
        coneIndex = connectedConeIDs(k);
        cP = squeeze(conePositions(coneIndex,:));
        coneSigma = coneSpacings(coneIndex)/3;
        coneProfile = exp(-0.5*((X-cP(1))/coneSigma).^2) .* exp(-0.5*((Y-cP(2))/coneSigma).^2);
        coneProfile(coneProfile>=flatTopZ) = flatTopZ;
        if (isempty(theRF))
            theRF = coneProfile * connectivityVectorForRGC(coneIndex);
        else
            theRF = theRF + coneProfile * connectivityVectorForRGC(coneIndex);
        end 
    end
    if (~isempty(theRF))
        theRF = theRF/max(theRF(:));
    end
end