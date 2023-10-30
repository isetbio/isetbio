function [theScenes, theNullStimulusScene, spatialSupportDegs] = ...
        generateStimulusFramesOnPresentationDisplay(...
                presentationDisplay, stimParams, ...
                spatialModulationPatterns, varargin)

    p = inputParser;
    p.addParameter('validateScenes', false, @islogical);
    p.addParameter('sceneIndexToCompute',  [], @isnumeric);
    p.parse(varargin{:});
    validateScenes = p.Results.validateScenes;
    sceneIndexToCompute = p.Results.sceneIndexToCompute;
    
    % Compute spatial support
    pixelsNum  = round(stimParams.stimSizeDegs / stimParams.pixelSizeDegs);
    spatialSupportDegs = linspace(-0.5*stimParams.stimSizeDegs, 0.5*stimParams.stimSizeDegs, pixelsNum);
    spatialSupportDegs = spatialSupportDegs - mean(spatialSupportDegs);

    % Background chromaticity and mean luminance vector
    xyY = [stimParams.backgroundChromaticity(1) stimParams.backgroundChromaticity(2) stimParams.backgroundLuminanceCdM2];
   
    % Background XYZ tri-stimulus values
    backgroundXYZ = (xyYToXYZ(xyY(:)))';
    
    % Background linear RGB primary values for the presentation display
    backgroundRGB = imageLinearTransform(backgroundXYZ, inv(displayGet(presentationDisplay, 'rgb2xyz')));
    
    % Background LMS excitations
    backgroundLMS = imageLinearTransform(backgroundRGB, displayGet(presentationDisplay, 'rgb2lms'));

    nStim = size(spatialModulationPatterns,1);

    theScenes = []; theNullStimulusScene = [];
    if (isempty(sceneIndexToCompute))
        theScenes = cell(1, nStim);
    end

    for sceneIndex = 0:nStim

        if (~isempty(sceneIndexToCompute))
            if (sceneIndex ~= sceneIndexToCompute)
                continue;
            end
        end

        % The LMS contrast image
        LMScontrastImage = zeros(size(spatialModulationPatterns,2), size(spatialModulationPatterns,2), numel(stimParams.coneContrasts));
        if (sceneIndex > 0)
            for coneIndex = 1:numel(stimParams.coneContrasts)
                LMScontrastImage(:,:,coneIndex) = stimParams.coneContrasts(coneIndex) * stimParams.contrast * squeeze(spatialModulationPatterns(sceneIndex,:,:));
            end
        end

        if (validateScenes)
            % Set the first pixel to 0 contrast so that we use it to estimate the
            % actual mean LMS excitation later on
            LMScontrastImage(1,1,:) = 0;
        end
        
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
        theScene = sceneFromFile(flipud(RGBsettings), format, meanLuminance, presentationDisplay);

        % Set the desired FOV 
        theScene = sceneSet(theScene, 'h fov', stimParams.stimSizeDegs);


        if (validateScenes)
            % Compute different scene representations for validation and visualization purposes
            sceneLMScontrastsImage = ...
                sceneRepresentations(theScene, presentationDisplay, backgroundLMS);

            % Assert that the scene cone contrasts match the desired ones
            figNo = 1999;
            assertDisplayContrasts(figNo, sceneLMScontrastsImage, LMScontrastImage);

            % Visualize scene components
            %figNo = 2000;
            %visualizeDisplayImage(figNo, sceneSRGBimage, sceneLMSexcitationsImage, presentationDisplay);
        end

        if (isempty(sceneIndexToCompute))
            if (sceneIndex == 0)
                theNullStimulusScene = theScene;
            else
                theScenes{sceneIndex} = theScene;
            end
        end

    end % for sceneIndex

    if (~isempty(sceneIndexToCompute))
         if (sceneIndexToCompute == 0)
             theNullStimulusScene = theScene;
         else
            theScenes{1} = theScene;
         end
    end

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

function assertDisplayContrasts(figNo, sceneLMScontrastsImage, desiredLMScontrastImage)
    hFig = figure(figNo); clf;
    targetRow = round(size(sceneLMScontrastsImage,1)/2);
    subplot(1,3,1)
    plot(sceneLMScontrastsImage(targetRow,:,1), 'r-'); hold on;
    plot(desiredLMScontrastImage(targetRow,:,1), 'k--');
    set(gca, 'YLim', [-0.15 0.15]);
    legend({'achieved', 'target'});
    title('L-cone contrast');

    subplot(1,3,2)
    plot(sceneLMScontrastsImage(targetRow,:,2), 'g-'); hold on;
    plot(desiredLMScontrastImage(targetRow,:,2), 'k--');
    set(gca, 'YLim', [-0.15 0.15]);
    legend({'achieved', 'target'});
    title('M-cone contrast');

    subplot(1,3,3)
    plot(sceneLMScontrastsImage(targetRow,:,3), 'b-'); hold on;
    plot(desiredLMScontrastImage(targetRow,:,3), 'k--');
    maxScontrast = max(squeeze(abs(desiredLMScontrastImage(targetRow,:,3))));
    if (maxScontrast < 0.15)
        maxScontrast = 0.15;
    end

    set(gca, 'YLim', [-1 1]*maxScontrast);
    legend({'achieved', 'target'});
    title('S-cone contrast');
end

function sceneLMScontrastsImage  = ...
    sceneRepresentations(theScene, presentationDisplay, backgroundLMS)

    emittedRadianceImage = sceneGet(theScene, 'energy');
    displayWavelengths = displayGet(presentationDisplay, 'wave');

    % Load the 2-deg Stockman cone fundamentals on a wavelength support matching the display
    coneFundamentals = ieReadSpectra(fullfile(isetbioDataPath,'human','stockman'), displayWavelengths);

    % Compute the LMS cone contrasts of the emitted radiance image
    sceneLMScontrastsImage = computeLMScontrastImage(emittedRadianceImage, coneFundamentals, backgroundLMS);
end

function LMScontrastImage = computeLMScontrastImage(radianceImage, coneFundamentals, coneExcitationsBackground)
    rowsNum = size(radianceImage,1);
    colsNum = size(radianceImage,2);
    wavelengthsNum = size(radianceImage, 3);
    radianceImage = reshape(radianceImage, [rowsNum*colsNum wavelengthsNum]);
        
    coneExcitationsImage = radianceImage * coneFundamentals;
    if (isempty(coneExcitationsBackground))
        % Assuming the top-left pixel has 0 cone contrast, i.e. equal to the background
        coneExcitationsBackground = coneExcitationsImage(1,:);
    end
    
    LMScontrastImage = bsxfun(@times, ...
        bsxfun(@minus, coneExcitationsImage, coneExcitationsBackground), ...
        1./coneExcitationsBackground);
    
    LMScontrastImage = reshape(LMScontrastImage, [rowsNum colsNum 3]);
end
