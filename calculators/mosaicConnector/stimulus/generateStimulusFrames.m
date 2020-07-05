function [theSceneFrames, presentationDisplay, sceneLuminanceSlice] = generateStimulusFrames(stimColor, stimSpatialParams, wavelengthSampling)

    % Generate the presentation display
    presentationDisplay = generatePresentationDisplay('wave', wavelengthSampling);

    % Background XYZ tri-stimulus values
    background.xyY = [stimColor.backgroundChroma(1) stimColor.backgroundChroma(2) stimColor.meanLuminanceCdPerM2];
   
    % Background XYZ tri-stimulus values
    background.XYZ = (xyYToXYZ(background.xyY(:)))';
    
    % Background linear RGB primary values for the presentation display
    background.RGB = imageLinearTransform(background.XYZ, inv(displayGet(presentationDisplay, 'rgb2xyz')));
    
    % Background LMS excitations
    background.LMS = imageLinearTransform(background.RGB, displayGet(presentationDisplay, 'rgb2lms'));
    
    % Generate a new scene for each spatial phase
    spatialPhases = (0 : stimSpatialParams.deltaPhaseDegs : (360-stimSpatialParams.deltaPhaseDegs));
    spatialPhases = spatialPhases/180*pi;
    theSceneFrames = cell(1, numel(spatialPhases));
    
    for spatialPhaseIndex = 1:numel(spatialPhases)
        % Stimulus spatial modulation of the L-, M-, and S-cone contrast
        test.LMScontrastImage = generateSpatialContrastImage(stimSpatialParams, stimColor.lmsContrast, spatialPhases(spatialPhaseIndex));
        
        % Stimulus LMS excitations image for the given background and spatial modulation
        test.LMSexcitationImage = bsxfun(@times, (1+test.LMScontrastImage), reshape(background.LMS, [1 1 3]));
    
        % Stimulus linear RGB primaries image
        test.RGBimage = imageLinearTransform(test.LMSexcitationImage, inv(displayGet(presentationDisplay, 'rgb2lms')));
        
        % Make sure we are in gamut (no subpixels with primary values outside of [0 1]
        outOfGamutPixels = numel(find((test.RGBimage(:)<0)|(test.RGBimage(:)>1)));
        assert(outOfGamutPixels==0, ...
            sprintf('%d subpixels with primary values > 1; %d subpixels with primary values < 0', ...
            numel(find(test.RGBimage>1)), numel(find(test.RGBimage<0))));
        
        % Generate a gamma corrected RGB image that we can pop in the
        % isetbio scene straightforward
        test.RGBimageGammaCorrected = (ieLUTLinear(test.RGBimage, displayGet(presentationDisplay, 'inverse gamma'))) / displayGet(presentationDisplay, 'nLevels');
    
        % Generate scene corresponding to the test stimulus on the presentation display
        format = 'rgb';
        meanLuminance = []; % DO NOT SET stimColor.meanLuminanceCdPerM2;
        theScene = sceneFromFile(test.RGBimageGammaCorrected, format, meanLuminance, presentationDisplay);
        
        % Set the desired FOV
        theScene = sceneSet(theScene, 'h fov', stimSpatialParams.fovDegs);
        
        % Add to the list of frames
        theSceneFrames{spatialPhaseIndex} = theScene;
           
        % Get the scene luminance
        luminanceImage = sceneGet(theScene, 'luminance');
        sceneLuminanceSlice(spatialPhaseIndex,:) = squeeze(luminanceImage(round(size(luminanceImage,1)/2),:));
        
        validateScene = ~true;
        if (validateScene)
            % Compute different scene representations for validation and visualization purposes
            [sceneLRGBimage, sceneSRGBimage, sceneLMScontrastsImage, sceneLMSexcitationsImage] = ...
                sceneRepresentations(theScene, presentationDisplay);

            % Assert that the scene cone contrasts match the desired ones
            figNo = 1999;
            assertDisplayContrasts(figNo, sceneLMScontrastsImage, test.LMScontrastImage);

            % Visualize scene components
            figNo = 2000;
            visualizeDisplayImage(figNo, sceneSRGBimage, sceneLMSexcitationsImage, presentationDisplay);
        end
        
    end
   
    
end


function LMScontrastImage = generateSpatialContrastImage(stimSpatialParams, lmsContrast, spatialPhase)
    x = linspace(-stimSpatialParams.fovDegs/2, stimSpatialParams.fovDegs/2, stimSpatialParams.pixelsNum);
    [X,Y] = meshgrid(x);
    
    fx = stimSpatialParams.gaborSpatialFrequencyCPD * cosd(stimSpatialParams.gaborOrientationDegs);
    fy = stimSpatialParams.gaborSpatialFrequencyCPD * sind(stimSpatialParams.gaborOrientationDegs);
    contrastPattern = cos(2*pi*( fx*(X-stimSpatialParams.gaborPosDegs(1)) + fy*(Y-stimSpatialParams.gaborPosDegs(2))) + spatialPhase);
    
    if (isinf(stimSpatialParams.gaborSigmaDegs))
        contrastEnvelope = ones(size(X));
    else
        if (numel(stimSpatialParams.gaborSigmaDegs) == 1)
            stimSpatialParams.gaborSigmaDegs = stimSpatialParams.gaborSigmaDegs*[1 1];
        end
        xx = X-stimSpatialParams.gaborPosDegs(1);
        yy = Y-stimSpatialParams.gaborPosDegs(2);
        contrastEnvelope = exp(-(0.5*(xx/stimSpatialParams.gaborSigmaDegs(1)).^2)) .* ...
                           exp(-(0.5*(yy/stimSpatialParams.gaborSigmaDegs(2)).^2));
    end
    
    LMScontrastImage = zeros(size(X,1), size(X,2), numel(lmsContrast));
    for coneIndex = 1:numel(lmsContrast) 
        LMScontrastImage(:,:,coneIndex) = lmsContrast(coneIndex) * contrastPattern .* contrastEnvelope;
    end
    
    % Set the first pixel to 0 contrast so that we use it to estimate the
    % actual mean LMS excitation later on
    LMScontrastImage(1,1,:) = 0;
end

function assertDisplayContrasts(figNo, sceneLMScontrastsImage, desiredLMScontrastImage)
    hFig = figure(figNo); clf;
    subplot(1,3,1)
    plot(sceneLMScontrastsImage(256,:,1), 'r-'); hold on;
    plot(desiredLMScontrastImage(256,:,1), 'k--');
    set(gca, 'YLim', [-0.15 0.15]);
    subplot(1,3,2)
    plot(sceneLMScontrastsImage(256,:,2), 'g-'); hold on;
    plot(desiredLMScontrastImage(256,:,2), 'k--');
    set(gca, 'YLim', [-0.15 0.15]);
    subplot(1,3,3)
    plot(sceneLMScontrastsImage(256,:,3), 'b-'); hold on;
    plot(desiredLMScontrastImage(256,:,3), 'k--');
    set(gca, 'YLim', [-0.15 0.15]);
end

function visualizeDisplayImage(figNo, sRGBimage, LMSimage,  presentationDisplay)
    backgroundLMS = LMSimage(1,1,:);
    %Compute the RGB image of the L-excitations image component
    tmp = LMSimage;
    for k = [2 3]
        tmp(:,:,k) = 0*tmp(:,:,k)+backgroundLMS(k);
    end
    tmp = imageLinearTransform(tmp, inv(displayGet(presentationDisplay, 'rgb2lms')));
    if ((min(tmp(:))<0) || (max(tmp(:))>1))
        tmp(tmp<0) = 0;
        tmp(tmp>1) = 1;
        fprintf('L-only component not realizable. Will clip.\n');
    end
    LexcitationSRGBimage = lrgb2srgb(tmp);
    
    % Compute the RGB image of the M-excitations image component
    tmp = LMSimage;
    for k = [1 3]
        tmp(:,:,k) = 0*tmp(:,:,k)+backgroundLMS(k);
    end
    tmp = imageLinearTransform(tmp, inv(displayGet(presentationDisplay, 'rgb2lms')));
    if ((min(tmp(:))<0) || (max(tmp(:))>1))
        tmp(tmp<0) = 0;
        tmp(tmp>1) = 1;
        fprintf('M-only component not realizable. Will clip.\n');
    end
    MexcitationSRGBimage = lrgb2srgb(tmp);
    

    % Compute the RGB image of the S-excitations image component
    tmp = LMSimage;
    for k =[1 2]
        tmp(:,:,k) = 0*tmp(:,:,k)+backgroundLMS(k);
    end
    tmp = imageLinearTransform(tmp, inv(displayGet(presentationDisplay, 'rgb2lms')));
    if ((min(tmp(:))<0) || (max(tmp(:))>1))
        tmp(tmp<0) = 0;
        tmp(tmp>1) = 1;
        fprintf('S-only component not realizable. Will clip.\n');
    end
    SexcitationSRGBimage = lrgb2srgb(tmp);
    

    hFig = figure(figNo); clf;

    % Plot the L-cone component
    theCurrentAxes = subplot(2,2,1); 
    image(theCurrentAxes, LexcitationSRGBimage);
    axis(theCurrentAxes, 'square');
    set(theCurrentAxes, 'XTick', [], 'YTick', []);
    title(theCurrentAxes, 'L-cone stimulus component');
    
    % Plot the M-cone component
    theCurrentAxes = subplot(2,2,2);
    image(theCurrentAxes, MexcitationSRGBimage);
    axis(theCurrentAxes, 'square');
    set(theCurrentAxes, 'XTick', [], 'YTick', []);
    title(theCurrentAxes, 'M-cone stimulus component');
     
    % Plot the S-cone component
    theCurrentAxes = subplot(2,2,3);
    image(theCurrentAxes, SexcitationSRGBimage);
    axis(theCurrentAxes, 'square');
    set(theCurrentAxes, 'XTick', [], 'YTick', []);
    title(theCurrentAxes, 'S-cone stimulus component');
    
    % Plot the composite stimulus
    theCurrentAxes = subplot(2,2,4);
    image(theCurrentAxes, sRGBimage);
    axis(theCurrentAxes, 'square')
    set(theCurrentAxes, 'XTick', [], 'YTick', []);
    title(theCurrentAxes, 'composite stimulus');
    drawnow
end




function [sceneLRGBimage, sceneSRGBimage, sceneLMScontrastsImage, sceneLMSexcitationsImage] = ...
    sceneRepresentations(theScene, presentationDisplay)

    emittedRadianceImage = sceneGet(theScene, 'energy');
    displaySPDs = displayGet(presentationDisplay, 'spd');
    displayWavelengths = displayGet(presentationDisplay, 'wave');
    [sceneLRGBimage, sceneSRGBimage] = displayRadianceToDisplayRGB(emittedRadianceImage, displaySPDs);
        
    % Load the 2-deg Stockman cone fundamentals on a wavelength support matching the display
    coneFundamentals = ieReadSpectra(fullfile(isetbioDataPath,'human','stockman'), displayWavelengths);

    % Compute the LMS cone contrasts of the emitted radiance image
    [sceneLMScontrastsImage, sceneLMSexcitationsImage] = displayRadianceToLMS(emittedRadianceImage, coneFundamentals);
end
    
function [LMScontrastsImage, LMSexcitationsImage] = displayRadianceToLMS(radianceImage, coneFundamentals)
    rowsNum = size(radianceImage,1);
    colsNum = size(radianceImage,2);
    wavelengthsNum = size(radianceImage, 3);
    radianceImage = reshape(radianceImage, [rowsNum*colsNum wavelengthsNum]);
        
    LMSexcitations = radianceImage * coneFundamentals;
    backgroundLMS = LMSexcitations(1,:);
    
    LMScontrasts = bsxfun(@times, ...
        bsxfun(@minus, LMSexcitations, backgroundLMS), ...
        1./backgroundLMS);
    
    LMScontrastsImage = reshape(LMScontrasts, [rowsNum colsNum 3]);
    LMSexcitationsImage = reshape(LMSexcitations, [rowsNum colsNum 3]);
end


