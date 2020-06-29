function theSceneFrames = generateStimulusFrames(stimColor, stimSpatialParams, wavelengthSampling)

    verifyAchievedContrasts = false;
    
    % Presentation display
    presentationDisplay = generatePresentationDisplay(...
        'viewingDistance', 0.57, ...
        'bitDepth', 8, ...
        'wave', wavelengthSampling);
    
    % Background XYZ tri-stimulus values
    background.XYZ = (xyYToXYZ([stimColor.backgroundChroma(:); stimColor.meanLuminanceCdPerM2]))';
    
    % Background linear RGB primary values for the presentation display
    background.RGB = imageLinearTransform(background.XYZ, inv(displayGet(presentationDisplay, 'rgb2xyz')));
    
    % Background LMS excitations
    background.LMS = imageLinearTransform(background.RGB, displayGet(presentationDisplay, 'rgb2lms'));
    
    % Generate a new scene for each spatial phase
    spatialPhases = (0 : stimSpatialParams.deltaPhaseDegs : (360-stimSpatialParams.deltaPhaseDegs))/180*pi;
    theSceneFrames = cell(1, numel(spatialPhases));
    
    for spatialPhaseIndex = 1:numel(spatialPhases)
        % Stimulus spatial modulation of the L-, M-, and S-cone contrast
        test.LMScontrastImage = generateSpatialContrastPatterns(stimSpatialParams, stimColor.lmsContrast, spatialPhases(spatialPhaseIndex));
        
        % Stimulus LMS excitations image for the given background and spatial modulation
        test.LMSexcitationImage = bsxfun(@times, (1+test.LMScontrastImage), reshape(background.LMS, [1 1 3]));
    
        % Stimulus linear RGB primaries image
        test.RGBimage = imageLinearTransform(test.LMSexcitationImage, inv(displayGet(presentationDisplay, 'rgb2lms')));

        % Make sure we are in gamut (no subpixels with primary values outside of [0 1]
        assert((numel(find(test.RGBimage>1))==0)&&(numel(find(test.RGBimage<0))==0), ...
            sprintf('%d subpixels with primary values > 1; %d subpixels with primary values < 0', ...
            numel(find(test.RGBimage>1)), numel(find(test.RGBimage<0))));
    
        % Gamma correct linear RGB values (primaries) through the display's
        % gamma table to get RGB settings values
        test.RGBimageGammaCorrected = ieLUTLinear(test.RGBimage, ...
            displayGet(presentationDisplay, 'inverse gamma'));
        test.RGBimageGammaCorrected = test.RGBimageGammaCorrected / max(test.RGBimageGammaCorrected(:));

        % Generate scene corresponding to the test stimulus on the presentation display
        theScene = sceneFromFile(test.RGBimageGammaCorrected,'rgb',stimColor.meanLuminanceCdPerM2, presentationDisplay);
        theScene = sceneSet(theScene, 'h fov', stimSpatialParams.fovDegs);

        % Add to the sceneFrames
        theSceneFrames{spatialPhaseIndex} = theScene;
        
        if (verifyAchievedContrasts)
            % Verify that we have the desired cone contrast profiles
            % Get the emitted radiance image
            emittedRadianceImage = sceneGet(theScene, 'energy');

            % Load the 2-deg Stockman cone fundamentals on a wavelength support matching the display
            coneFundamentals = ieReadSpectra(fullfile(isetbioDataPath,'human','stockman'), displayGet(presentationDisplay, 'wave'));

            % Compute the LMS cone contrasts of the emitted radiance image
            test.achievedLMScontrastImage = computeLMScontrastImage(emittedRadianceImage, coneFundamentals, background.LMS);
            
            diffs = test.LMScontrastImage - test.achievedLMScontrastImage;
            maxLconeContrast = max(max(squeeze(test.LMScontrastImage(:,:,1))));
            maxMconeContrast = max(max(squeeze(test.LMScontrastImage(:,:,2))));
            maxSconeContrast = max(max(squeeze(test.LMScontrastImage(:,:,3))));

            maxLconeContrastDeviation = max(max(abs(squeeze(diffs(:,:,1)))));
            maxMconeContrastDeviation = max(max(abs(squeeze(diffs(:,:,2)))));
            maxSconeContrastDeviation = max(max(abs(squeeze(diffs(:,:,3)))));

            fprintf('Contrast deviations: %2.3f%% (L), %2.3f%% (M), %2.3f%% (S)\n', ...
                maxLconeContrastDeviation/maxLconeContrast*100, ...
                maxMconeContrastDeviation/maxMconeContrast*100, ...
                maxSconeContrastDeviation/maxSconeContrast*100 ...
                );
        end
        
    end
    
end

function LMScontrastImage = generateSpatialContrastPatterns(stimSpatialParams, lmsContrast, spatialPhase)
    x = linspace(-stimSpatialParams.fovDegs/2, stimSpatialParams.fovDegs/2, stimSpatialParams.pixelsNum);
    [X,Y] = meshgrid(x);
    
    fx = stimSpatialParams.gaborSpatialFrequencyCPD * cosd(stimSpatialParams.gaborOrientationDegs);
    fy = stimSpatialParams.gaborSpatialFrequencyCPD * sind(stimSpatialParams.gaborOrientationDegs);
    contrastPattern = cos(2*pi*( fx*(X-stimSpatialParams.gaborPosDegs(1)) + fy*(Y-stimSpatialParams.gaborPosDegs(2))) + spatialPhase);
    
    if (isinf(stimSpatialParams.gaborSigma))
        contrastEnvelope = ones(size(X));
    else
        if (numel(stimSpatialParams.gaborSigma) == 1)
            stimSpatialParams.gaborSigma = stimSpatialParams.gaborSigma*[1 1];
        end
        contrastEnvelope = exp(-((X-stimSpatialParams.gaborPosDegs(1))/stimSpatialParams.gaborSigma(1)).^2) .* exp(-((Y-stimSpatialParams.gaborPosDegs(2))/stimSpatialParams.gaborSigma(2)).^2);
    end
    
    LMScontrastImage = zeros(size(X,1), size(X,2), numel(lmsContrast));
    for coneIndex = 1:numel(lmsContrast) 
        LMScontrastImage(:,:,coneIndex) = lmsContrast(coneIndex) * contrastPattern .* contrastEnvelope;
    end
end

function LMScontrastImage = computeLMScontrastImage(radianceImage, coneFundamentals, coneExcitationsBackground)
    rowsNum = size(radianceImage,1);
    colsNum = size(radianceImage,2);
    wavelengthsNum = size(radianceImage, 3);
    radianceImage = reshape(radianceImage, [rowsNum*colsNum wavelengthsNum]);
        
    coneExcitationsImage = radianceImage * coneFundamentals;

    LMScontrastImage = bsxfun(@times, ...
        bsxfun(@minus, coneExcitationsImage, coneExcitationsBackground), ...
        1./coneExcitationsBackground);
    
    LMScontrastImage = reshape(LMScontrastImage, [rowsNum colsNum 3]);
end

