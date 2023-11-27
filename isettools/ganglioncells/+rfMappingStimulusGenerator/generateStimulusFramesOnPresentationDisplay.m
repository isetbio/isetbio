function [theScenes, theNullStimulusScene, spatialSupportDegs, coneFundamentalsStruct] = ...
        generateStimulusFramesOnPresentationDisplay(...
                presentationDisplay, stimParams, ...
                spatialModulationPatterns, varargin)

    p = inputParser;
    p.addParameter('validateScenes', false, @islogical);
    p.addParameter('sceneIndexToCompute',  [], @isnumeric);
    p.addParameter('customConeFundamentals', [], @(x)(isempty(x)||isstruct(x)));
    p.addParameter('withPreviouslyComputedConeFundamentalsStruct', [], @(x)(isempty(x)||isstruct(x)));
    p.parse(varargin{:});
    validateScenes = p.Results.validateScenes;
    sceneIndexToCompute = p.Results.sceneIndexToCompute;
    customConeFundamentals = p.Results.customConeFundamentals;
    previouslyComputedConeFundamentalsStruct = p.Results.withPreviouslyComputedConeFundamentalsStruct;

    if (isempty(previouslyComputedConeFundamentalsStruct))
        if (isempty(customConeFundamentals))
            % Load the 2-deg Stockman cone fundamentals on wavelength support matching the display
            displayWavelengths = displayGet(presentationDisplay, 'wave');
            coneFundamentals = ieReadSpectra(fullfile(isetbioDataPath,'human','stockman'), displayWavelengths);

            % Returned coneFundamentals struct
            coneFundamentalsStruct.coneFundamentals = coneFundamentals;
            coneFundamentalsStruct.spectralSupport = displayWavelengths;

            fprintf('Employing the standard SS-2 cone fundamentals.\n');
        else
            displayWavelengths = displayGet(presentationDisplay, 'wave');
            assert(isfield(customConeFundamentals, 'wavelengthSupport'), ...
                'customConeFundamentals does not contain wavelength support info');
            assert(isfield(customConeFundamentals, 'quantalExcitationSpectra'), ...
                'customConeFundamentals does not contain quantalExcitationSpectra info');
            assert(size(customConeFundamentals.quantalExcitationSpectra,2) == 3, ...
                'customConeFundamentals.spd is not an Nx3 matrix');
            assert(size(customConeFundamentals.quantalExcitationSpectra,1) == numel(customConeFundamentals.wavelengthSupport), ...
                'customConeFundamentals.spf does not have the same dimensionality as customConeFundamentals.wavelengthSupport');
    
            if (~isequal(displayWavelengths, customConeFundamentals.wavelengthSupport))
                % Resample customConeFundamentals.spd to wavelength support matching the display
                resampledCustomConeFundamantals = displayWavelengths*0;
                for iChannel = 1:size(customConeFundamentals.quantalExcitationSpectra,2)
                    resampledCustomConeFundamantals(:,iChannel) = interp1(...
                        customConeFundamentals.wavelengthSupport, customConeFundamentals.quantalExcitationSpectra(:,iChannel), ...
                        displayWavelengths, 'linear','extrap');
                end
                customConeFundamentals.quantalExcitationSpectra = resampledCustomConeFundamantals;
                customConeFundamentals.wavelengthSupport = displayWavelengths;
            end
            
            coneFundamentals = customConeFundamentals.quantalExcitationSpectra/max(customConeFundamentals.quantalExcitationSpectra(:));

            % Returned coneFundamentals struct
            coneFundamentalsStruct.coneFundamentals = coneFundamentals;
            coneFundamentalsStruct.spectralSupport = displayWavelengths;

            % Compare to default SS2
            compareCustomConeFundamentalsToDefaultSS2(customConeFundamentals);

            fprintf('Employing the custom cone fundamentals.\n');
        end
    else
        coneFundamentals = previouslyComputedConeFundamentalsStruct.coneFundamentals;
        displayWavelengths = previouslyComputedConeFundamentalsStruct.spectralSupport;
        coneFundamentalsStruct = previouslyComputedConeFundamentalsStruct;
        fprintf('Employing previously used cone fundamentals.\n');
    end

    % Compute the displayRGCtoLMS matrix
    displayRGBtoLMS = (coneFundamentals' * displayGet(presentationDisplay, 'spd', displayWavelengths))';
    displayLMStoRGB = inv(displayRGBtoLMS);

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
    backgroundLMS = imageLinearTransform(backgroundRGB, displayRGBtoLMS);

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
        RGBimage = imageLinearTransform(LMSexcitationImage, displayLMStoRGB);

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
            emittedRadianceImage = sceneGet(theScene, 'energy');
            % Compute the LMS cone contrasts of the emitted radiance image
            sceneLMScontrastsImage = computeLMScontrastImage(emittedRadianceImage, coneFundamentals, backgroundLMS);
       

            % Assert that the scene cone contrasts match the desired ones
            figNo = 1999;
            assertDisplayContrasts(figNo, sceneLMScontrastsImage, LMScontrastImage);

            % Visualize scene components
            figure(2000);
            ax = subplot(1,1,1);
            visualizeScene(theScene, ...
                'presentationDisplay', presentationDisplay, ...
                'axesHandle', ax);
            drawnow;
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


function compareCustomConeFundamentalsToDefaultSS2(customConeFundamentals)
    % Normalize
    coneFundamentals = customConeFundamentals.quantalExcitationSpectra/max(customConeFundamentals.quantalExcitationSpectra(:));
    StockmanSharpe2DegConeFundamentals = ieReadSpectra(fullfile(isetbioDataPath,'human','stockman'), customConeFundamentals.wavelengthSupport);

    hFig = figure(223); clf;
    subplot(2,2,1);
    plot(customConeFundamentals.wavelengthSupport, StockmanSharpe2DegConeFundamentals(:,1), 'r--', 'LineWidth', 1.5);
    hold on;
    plot(customConeFundamentals.wavelengthSupport, StockmanSharpe2DegConeFundamentals(:,2), 'g--', 'LineWidth', 1.5);
    plot(customConeFundamentals.wavelengthSupport, StockmanSharpe2DegConeFundamentals(:,3), 'b--', 'LineWidth', 1.5);
    title('Stockman 2 deg cone fundamentals')

    subplot(2,2,2);
    plot(customConeFundamentals.wavelengthSupport, customConeFundamentals.quantalExcitationSpectra(:,1), 'r-', 'LineWidth', 1.5);
    hold on;
    plot(customConeFundamentals.wavelengthSupport, customConeFundamentals.quantalExcitationSpectra(:,2), 'g-', 'LineWidth', 1.5);
    plot(customConeFundamentals.wavelengthSupport, customConeFundamentals.quantalExcitationSpectra(:,3), 'b-', 'LineWidth', 1.5);
    title('cMosaic cone fundamentals')

    subplot(2,2,3);
    plot(customConeFundamentals.wavelengthSupport, customConeFundamentals.quantalExcitationSpectra(:,1), 'r-', 'LineWidth', 1.5);
    hold on;
    plot(customConeFundamentals.wavelengthSupport, customConeFundamentals.quantalExcitationSpectra(:,2), 'g-', 'LineWidth', 1.5);
    plot(customConeFundamentals.wavelengthSupport, customConeFundamentals.quantalExcitationSpectra(:,3), 'b-', 'LineWidth', 1.5);
    plot(customConeFundamentals.wavelengthSupport, StockmanSharpe2DegConeFundamentals(:,1), 'k--', 'LineWidth', 1.5);
    plot(customConeFundamentals.wavelengthSupport, StockmanSharpe2DegConeFundamentals(:,2), 'k--', 'LineWidth', 1.5);
    plot(customConeFundamentals.wavelengthSupport, StockmanSharpe2DegConeFundamentals(:,3), 'k--', 'LineWidth', 1.5);
    title('cMosaic cone fundamentals')

    subplot(2,2,4);
    plot(customConeFundamentals.wavelengthSupport, customConeFundamentals.quantalExcitationSpectra(:,1)./StockmanSharpe2DegConeFundamentals(:,1), 'r-', 'LineWidth', 1.5);
    hold on;
    plot(customConeFundamentals.wavelengthSupport, customConeFundamentals.quantalExcitationSpectra(:,2)./StockmanSharpe2DegConeFundamentals(:,2), 'g-', 'LineWidth', 1.5);
    plot(customConeFundamentals.wavelengthSupport, customConeFundamentals.quantalExcitationSpectra(:,3)./StockmanSharpe2DegConeFundamentals(:,3), 'b-', 'LineWidth', 1.5);
    plot(customConeFundamentals.wavelengthSupport, customConeFundamentals.wavelengthSupport*0 + 1, 'k-');
    title('cMosaic cone fundamentals ./ SS2')
    drawnow;
end
