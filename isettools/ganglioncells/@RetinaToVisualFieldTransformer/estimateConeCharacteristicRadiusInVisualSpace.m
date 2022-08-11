function dStruct = estimateConeCharacteristicRadiusInVisualSpace(obj,...
            whichEye, eccDegs, testSubjectID, pupilDiameterMM,  dataFileName, varargin)
% Estimate the characteristic radius of the mean cone at the center of a
% cone mosaic centered at a given (x,y) eccentricity in degrees
% as projected on to visual space using optics at the same eccentricity

    % Parse input
    p = inputParser;
    p.addParameter('anatomicalConeCharacteristicRadiusDegs', [], @(x)(isempty(x)||isscalar(x)));
    p.addParameter('hFig', []);
    p.addParameter('videoOBJ', []);
    p.parse(varargin{:});

    anatomicalConeCharacteristicRadiusDegs = p.Results.anatomicalConeCharacteristicRadiusDegs;
    hFig = p.Results.hFig;
    videoOBJ = p.Results.videoOBJ;

    % Cone aperture modifiers
    % The cone light gathering aperture is described by a Gaussian function with 
    % sigma equal to 0.204 x inner segment diameter (cone diameter)
    sigmaGaussian = 0.204;  % From McMahon et al, 2000

    % Set cone aperture modifiers
    coneApertureModifiers = struct(...
        'smoothLocalVariations', true, ...
        'sigma',  sigmaGaussian, ...
        'shape', 'Gaussian');

    % Generate the cone mosaic at the given eye and eccentricity
    cm = cMosaic(...
        'whichEye', whichEye, ...
        'sizeDegs', [1 1], ...
        'eccentricityDegs', eccDegs, ...
        'rodIntrusionAdjustedConeAperture', true, ...
        'customDegsToMMsConversionFunction', @RGCmodels.Watson.convert.rhoDegsToMMs, ...
        'customMMsToDegsConversionFunction', @RGCmodels.Watson.convert.rhoMMsToDegs, ...
        'wave', 550);

    % How many cones we have in this [1x1] deg patch
    conesNumInRetinalPatch = numel(cm.coneTypes);

    % If an anatomicalConeCharacteristicRadiusDegs) is not specified,
    % compute it here from the center of the mosaic
    if (isempty(anatomicalConeCharacteristicRadiusDegs))
        % If we have less than 6 cones, we will not do any analysis
        if (conesNumInRetinalPatch < 6)
            fprintf(2, 'Less than 6 cones at this eccentricity, skipping computation of cone aperture in visual space.\n')
            
            % Return struct
            dStruct.conesNumInRetinalPatch = conesNumInRetinalPatch;
            dStruct.anatomicalConeCharacteristicRadiusDegs = nan;
            dStruct.visualConeCharacteristicRadiusDegs = nan;
            return;
        end

        % Sort cones according to their distance from the mosaic center
        coneDistancesFromMosaicCenter = sqrt(sum(bsxfun(@minus, cm.coneRFpositionsDegs, cm.eccentricityDegs).^2,2));
        [~,idx] = sort(coneDistancesFromMosaicCenter, 'ascend');

        % Estimate mean anatomical cone aperture in the mosaic'c center
        meanConeApertureDegsInMosaicCenter = mean(cm.coneApertureDiametersDegs(idx(1:6)));
        anatomicalConeCharacteristicRadiusDegs = 0.204 * sqrt(2.0) * meanConeApertureDegsInMosaicCenter;
    end


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
        % Return struct
        dStruct.conesNumInRetinalPatch = conesNumInRetinalPatch;
        dStruct.anatomicalConeCharacteristicRadiusDegs = anatomicalConeCharacteristicRadiusDegs;
        dStruct.visualConeCharacteristicRadiusDegs = nan;
        return;
    end

    % Extract the PSF
    thePSFData = psfEnsemble{1};
    targetWavelength = cm.wave(1);
    [~, wIndex] = min(abs(thePSFData.supportWavelength-targetWavelength));
    thePSFData.data = squeeze(thePSFData.data(:,:,wIndex));


    if (~isempty(dataFileName))
        p = getpref('ISETMacaque');
        pdfFileName = fullfile(p.generatedDataDir, 'coneApertureBackingOut', ...
                sprintf('%s_ecc_%2.0f_%2.0f.pdf',dataFileName, cm.eccentricityDegs(1), cm.eccentricityDegs(2)));
    else
        pdfFileName = 'tmp';
    end

    % Compute the visually-projected cone aperture given this PSF
    visualConeCharacteristicRadiusDegs = RetinaToVisualFieldTransformer.analyzeEffectOfPSFonConeAperture(...
                        anatomicalConeCharacteristicRadiusDegs, thePSFData, ...
                        hFig, videoOBJ, pdfFileName);
    
    % Return struct
    dStruct.conesNumInRetinalPatch = conesNumInRetinalPatch;
    dStruct.anatomicalConeCharacteristicRadiusDegs = anatomicalConeCharacteristicRadiusDegs;
    dStruct.visualConeCharacteristicRadiusDegs = visualConeCharacteristicRadiusDegs;
end