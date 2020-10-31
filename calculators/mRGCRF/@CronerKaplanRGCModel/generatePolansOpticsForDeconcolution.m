function [theOptics, thePSF, thePSFsupportDegs] = generatePolansOpticsForDeconcolution(PolansSubjectID, imposedRefractionErrorDiopters, ...
    pupilDiameterMM , wavelengthSampling, micronsPerDegree, patchEcc, varargin)
    
     % Parse input
    p = inputParser;
    p.addParameter('eccentricityUnits', 'degrees', @(x)(ismember(x, {'microns', 'degrees'})));
    p.addParameter('noLCA', ~true, @islogical);
    p.addParameter('noOptics', ~true, @islogical);
    p.addParameter('doNotZeroCenterPSF', false, @islogical);
    p.parse(varargin{:});
    
    eccentricityUnits = p.Results.eccentricityUnits;
    noLCA = p.Results.noLCA;
    noOptics = p.Results.noOptics;
    doNotZeroCenterPSF = p.Results.doNotZeroCenterPSF;
    
    if (strcmp(eccentricityUnits, 'degrees'))
        patchEccDegs = patchEcc;
    else
        patchEccDegs = WatsonRGCModel.rhoMMsToDegs(patchEcc*1e-3);
    end
    
    % Choose wavefront sampled based on the eccentricity of the mosaic,
    % because PSFs get larger as ecc increases
    patchEccRadiusDegs = sqrt(sum(patchEccDegs.^2,2));
    
    if (patchEccRadiusDegs <= 10)
        wavefrontSpatialSamples = 501;
    elseif (patchEccRadiusDegs <= 15)
        wavefrontSpatialSamples = 801;
    else
        wavefrontSpatialSamples = 101;
    end
    wavefrontSpatialSamples = 1001;
    
    eccXrangeDegs = patchEccDegs(1)*[1 1];
    eccYrangeDegs = patchEccDegs(2)*[1 1];
    deltaEcc = 1;
    [~, ~, thePSFs, thePSFsupportDegs, theOIs] = CronerKaplanRGCModel.psfsAtEccentricity(PolansSubjectID, ...
                imposedRefractionErrorDiopters, pupilDiameterMM, wavelengthSampling, micronsPerDegree, ...
                wavefrontSpatialSamples, eccXrangeDegs, eccYrangeDegs, deltaEcc, ...
                'noLCA', noLCA, 'noOptics', noOptics, 'doNotZeroCenterPSF', doNotZeroCenterPSF);

    theOptics = theOIs{1,1,1};
    thePSF = thePSFs(1,1,1,:,:,:);

end
