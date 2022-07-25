function [theOTFData, thePSFData] = vLambdaWeightedPSFandOTF(obj, cm, testSubjectID, pupilDiameterMM)
    % Generate optics for this eye, eccentricity, subject, and pupil size
    switch (obj.ZernikeDataBase)
        % Artal
        case RetinaToVisualFieldTransformer.Artal
            subtractCentralRefraction = ArtalOptics.constants.subjectRequiresCentralRefractionCorrection(...
                cm.whichEye, testSubjectID);
        % Polans
        case RetinaToVisualFieldTransformer.Polans
            subtractCentralRefraction = PolansOptics.constants.subjectRequiresCentralRefractionCorrection(...
                cm.whichEye, testSubjectID);
    end

    [oiEnsemble, psfEnsemble] = cm.oiEnsembleGenerate(cm.eccentricityDegs, ...
                    'zernikeDataBase', obj.ZernikeDataBase, ...
                    'subjectID', testSubjectID, ...
                    'pupilDiameterMM', pupilDiameterMM, ...
                    'zeroCenterPSF', true, ...
                    'subtractCentralRefraction', subtractCentralRefraction, ...
                    'wavefrontSpatialSamples', 701, ...
                    'warningInsteadOfErrorForBadZernikeCoeffs', true);

    if (isempty(oiEnsemble))
        fprintf(2, 'Could not generate optics at this eccentricity.\n');
        % Return empties
        theOTFData = []; 
        thePSFData = [];
        return;
    end

    % Extract the OTF & the PSF
    thePSFData = psfEnsemble{1};
    theOI = oiEnsemble{1};
    theOTFData.data = theOI.optics.OTF.OTF;
    theOTFData.supportX = theOI.optics.OTF.fx;
    theOTFData.supportY = theOI.optics.OTF.fy;
    theOTFData.supportWavelength = theOI.optics.OTF.wave;

    % Compute v_lambda weights for weigthing the PSF/OTF
    weights = vLambdaWeigts(cm.wave);

    % Compute vLambda weighted OTF and PSF
    vLambdaWeightedOTF = zeros(size(theOTFData.data,1), size(theOTFData.data,2));
    vLambdaWeightedPSF = zeros(size(thePSFData.data,1), size(thePSFData.data,2));
    for iWave = 1:size(theOTFData.data,3)
        vLambdaWeightedOTF = vLambdaWeightedOTF + theOTFData.data(:,:,iWave) * weights(iWave);
        vLambdaWeightedPSF = vLambdaWeightedPSF + thePSFData.data(:,:,iWave) * weights(iWave);
    end
    theOTFData.data = vLambdaWeightedOTF;
    thePSFData.data = vLambdaWeightedPSF;

    % Remove support wavelength
    theOTFData = rmfield(theOTFData, 'supportWavelength');
    thePSFData = rmfield(thePSFData, 'supportWavelength');
    thePSFData = rmfield(thePSFData, 'zCoeffs');
end


function w = vLambdaWeigts(wavelengthSupport)
    load T_xyz1931;
    S = WlsToS(wavelengthSupport(:));
    T_vLambda = SplineCmf(S_xyz1931,T_xyz1931(2,:),S);
    w = T_vLambda/max(T_vLambda(:));
    w = w / sum(w(:));
end
