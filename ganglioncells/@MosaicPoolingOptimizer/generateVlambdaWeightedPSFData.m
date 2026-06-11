function thePSFData = generateVlambdaWeightedPSFData(theMidgetRGCMosaic, opticsParams)
    
    dataOut = theMidgetRGCMosaic.generateOptics(opticsParams);
    thePSFData = dataOut.thePSFData;

    % Load vLambda
    load T_xyz1931;
    S = WlsToS(thePSFData.supportWavelength);
    vLambda = SplineCmf(S_xyz1931, T_xyz1931(2,:), S);
    vLambda = vLambda / sum(vLambda);
    
    % Weight PSF by vLambda
    for iW = 1:numel(vLambda)
        if (iW == 1)
            thePSFmap = squeeze(thePSFData.data(:,:,iW)) * vLambda(iW);
        else
            thePSFmap = thePSFmap + squeeze(thePSFData.data(:,:,iW)) * vLambda(iW);
        end
    end

    thePSFData.vLambdaWeightedPSF = flipud(thePSFmap);
end