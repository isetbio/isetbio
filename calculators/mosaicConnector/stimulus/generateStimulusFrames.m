function [theSceneFrames, presentationDisplay, stimSpatialParamsPixelsNum] = generateStimulusFrames(stimColor, stimSpatialParams, wavelengthSampling)

    % Generate the presentation display
    presentationDisplay = generatePresentationDisplay('wave', wavelengthSampling);

    % Background chromaticity and mean luminance vector
    xyY = [stimColor.backgroundChroma(1) stimColor.backgroundChroma(2) stimColor.meanLuminanceCdPerM2];
   
    % Background XYZ tri-stimulus values
    backgroundXYZ = (xyYToXYZ(xyY(:)))';
    
    % Background linear RGB primary values for the presentation display
    backgroundRGB = imageLinearTransform(backgroundXYZ, inv(displayGet(presentationDisplay, 'rgb2xyz')));
    
    % Background LMS excitations
    backgroundLMS = imageLinearTransform(backgroundRGB, displayGet(presentationDisplay, 'rgb2lms'));
    
    % Generate a new scene for each spatial phase
    spatialPhases = (0 : stimSpatialParams.deltaPhaseDegs : (360-stimSpatialParams.deltaPhaseDegs));
    spatialPhases = spatialPhases/180*pi;
    theSceneFrames = cell(1, numel(spatialPhases));
    
    % Stimulus LMS contrast vector
    stimLMScontrast = stimColor.lmsContrast;
    
    % Compute pixelsNum
    pixelsNumAtCurrentSpatialFrequencyAndMinPixelsPerCycle = stimSpatialParams.gaborSpatialFrequencyCPD * stimSpatialParams.pixelsNumPerSF;
    
    % Store it in stimSpatialParams struct and also return it to the caller function
    fprintf('Pixels num: max([%2.0f %2.0f])\n', stimSpatialParams.minPixelsNum, pixelsNumAtCurrentSpatialFrequencyAndMinPixelsPerCycle);
    stimSpatialParams.pixelsNum = max([stimSpatialParams.minPixelsNum  pixelsNumAtCurrentSpatialFrequencyAndMinPixelsPerCycle]);
    stimSpatialParamsPixelsNum = stimSpatialParams.pixelsNum;

    % Generate scenes for each spatial phase
    parfor spatialPhaseIndex = 1:numel(spatialPhases)
        fprintf('Generating scene for the %2.0f deg spatial phase stimulus.\n', spatialPhases(spatialPhaseIndex)/pi*180);
        
        % Stimulus spatial modulation of the L-, M-, and S-cone contrast
        LMScontrastImage = generateSpatialContrastImage(stimSpatialParams, stimLMScontrast, spatialPhases(spatialPhaseIndex));
        
        % Stimulus LMS excitations image for the given background and spatial modulation
        LMSexcitationImage = bsxfun(@times, (1+LMScontrastImage), reshape(backgroundLMS, [1 1 3]));
    
        % Stimulus linear RGB primaries image
        RGBimage = imageLinearTransform(LMSexcitationImage, inv(displayGet(presentationDisplay, 'rgb2lms')));
        
        % Make sure we are in gamut (no subpixels with primary values outside of [0 1]
        outOfGamutPixels = numel(find((RGBimage(:)<0)|(RGBimage(:)>1)));
        assert(outOfGamutPixels==0, ...
            sprintf('%d subpixels with primary values > 1; %d subpixels with primary values < 0', ...
            numel(find(RGBimage>1)), numel(find(RGBimage<0))));
        
        % Generate a gamma corrected RGB image (RGBsettings) that we can pop in the
        % isetbio scene straightforward
        RGBsettings = (ieLUTLinear(RGBimage, displayGet(presentationDisplay, 'inverse gamma'))) / displayGet(presentationDisplay, 'nLevels');
    
        % Generate scene corresponding to the test stimulus on the presentation display
        format = 'rgb';
        meanLuminance = []; % EMPTY, so that mean luminance is determined from the rgb settings values we pass
        theScene = sceneFromFile(RGBsettings, format, meanLuminance, presentationDisplay);
        
        % Set the desired FOV
        theScene = sceneSet(theScene, 'h fov', stimSpatialParams.fovDegs);
        
        % Add to the list of frames
        theSceneFrames{spatialPhaseIndex} = theScene;
           
        validateScene = ~true;
        if (validateScene)
            % Compute different scene representations for validation and visualization purposes
            [sceneLRGBimage, sceneSRGBimage, sceneLMScontrastsImage, sceneLMSexcitationsImage] = ...
                sceneRepresentations(theScene, presentationDisplay);

            % Assert that the scene cone contrasts match the desired ones
            figNo = 1999;
            assertDisplayContrasts(figNo, sceneLMScontrastsImage, LMScontrastImage);

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
