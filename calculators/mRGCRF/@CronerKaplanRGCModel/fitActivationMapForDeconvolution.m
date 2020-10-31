function [rfSigmas, rfGain, hiResConeActivationMap, hiResPSFsupportDegs] = fitActivationMapForDeconvolution(functionName, coneActivations, conePosDegs, coneAperturesDegs, thePSFsupportDegs)

    deltaX = thePSFsupportDegs(2)-thePSFsupportDegs(1);
    minConeCharacteristicRadiusDegs = 0.5*min(coneAperturesDegs)/3;
    
    % Fit Gaussian to cone activations
    [fittedParams, rfFunction] = CronerKaplanRGCModel.fitElliptical2DGausianToRF(functionName, ...
        conePosDegs(:,1), conePosDegs(:,2), coneActivations(:), ...
        deltaX/40, minConeCharacteristicRadiusDegs, [0 0]);
    
    
    % Retrieve fit params
    rfGain = fittedParams(1);
    if (strcmp(functionName, 'elliptical Gaussian'))
        rfSigmas = [fittedParams(4) fittedParams(5)];
        fprintf('fitted elliptical gaussian exponent: %2.2f\n',  fittedParams(7));
    else
        rfSigmas = fittedParams(4)*[1 1];
        fprintf('fitted circular gaussian exponent: %2.2f\n',  fittedParams(5));
    end
    
    % Generate a high-res fitted map
    maxPSFsupport = deltaX*100;
    hiResPSFsupportDegs = -maxPSFsupport:deltaX:maxPSFsupport;
    [XXX,YYY] = meshgrid(hiResPSFsupportDegs, hiResPSFsupportDegs);
    xyData(:,1) = XXX(:);
    xyData(:,2) = YYY(:);
    fittedActivationMap = rfFunction(fittedParams,xyData);
    hiResConeActivationMap = reshape(fittedActivationMap, [numel(hiResPSFsupportDegs) numel(hiResPSFsupportDegs)]);
end
