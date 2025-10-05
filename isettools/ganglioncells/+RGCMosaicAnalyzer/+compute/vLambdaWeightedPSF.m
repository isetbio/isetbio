function theVlambdaWeightedPSF = vLambdaWeightedPSF(thePSF)
    % Load vLambda
    load T_xyz1931;
    S = WlsToS(thePSF.supportWavelength);
    vLambda = SplineCmf(S_xyz1931, T_xyz1931(2,:), S);
    vLambda = vLambda / sum(vLambda);
    
    % Weight PSF by vLambda
    for iW = 1:numel(vLambda)
        if (iW == 1)
            thePSFmap = squeeze(thePSF.data(:,:,iW)) * vLambda(iW);
        else
            thePSFmap = thePSFmap + squeeze(thePSF.data(:,:,iW)) * vLambda(iW);
        end
    end

    theVlambdaWeightedPSF = thePSFmap;
end

