function t_generateConeSpecificStimuli
% Illustrate how to generate scenes depicting various cone-specific stimuli
%
% Syntax:
%   t_generateConeSpecificStimuli
%
% Description:
%    Simple script that demonstrates how to generate scenes depicting 
%    various superimposed cone-specific stimuli on a specific background.
%
% Inputs:
%    None.
%
% Outputs:
%    None.
%
% Optional key/value pairs:
%    None.
%

% History:
%    03/28/20  npc  Wrote it.

    %% Specify the desired mean (x,y) chromaticity and luminance (background)
    % Horizon light - see https://en.wikipedia.org/wiki/Standard_illuminant
    background.xyChroma = [0.345 0.358];
    % 40 cd/m2
    background.luminance = 40;
    % In vector form
    background.xyY = [background.xyChroma(1) background.xyChroma(2) background.luminance];
    
    %% Specify stimulus type
    stimType = 'S-(L+M)';
    stimType = 'L-M';
    stimType = 'L+M';
    stimType = 'L+M+S';
    stimType = 'all-different';
    test.spatialModulationParams = getModulationParamsForStimType(stimType);
            
    %% Stimulus field of view
    fieldOfViewDegs = 4.0;
    
    %% Presentation display
    presentationDisplay = generatePresentationDisplay(); 
    
    %% Background XYZ tri-stimulus values
    background.XYZ = (xyYToXYZ(background.xyY(:)))';
    
    %% Background linear RGB primary values for the presentation display
    background.RGB = imageLinearTransform(background.XYZ, inv(displayGet(presentationDisplay, 'rgb2xyz')));
    
    %% Background LMS excitations
    background.LMS = imageLinearTransform(background.RGB, displayGet(presentationDisplay, 'rgb2lms'));
    
    %% Stimulus spatial modulation of the L-, M-, and S-cone contrast
    test.LMScontrastImage = generateSpatialContrastPatterns(fieldOfViewDegs, test.spatialModulationParams);
    
    %% Stimulus LMS excitations image for the given background and spatial modulation
    test.LMSexcitationImage = bsxfun(@times, (1+test.LMScontrastImage), reshape(background.LMS, [1 1 3]));
    
    %% Stimulus linear RGB primaries image
    test.RGBimage = imageLinearTransform(test.LMSexcitationImage, inv(displayGet(presentationDisplay, 'rgb2lms')));
    
    %% Make sure we are in gamut (no subpixels with primary values outside of [0 1]
    assert((numel(find(test.RGBimage>1))==0)&&(numel(find(test.RGBimage<0))==0), ...
        sprintf('%d subpixels with primary values > 1; %d subpixels with primary values < 0', ...
        numel(find(test.RGBimage>1)), numel(find(test.RGBimage<0))));
    
    %% Gamma correction
    % Gamma correct linear RGB values (primaries) through the display's
    % gamma table to get RGB settings values
    test.RGBimageGammaCorrected = ieLUTLinear(test.RGBimage, ...
        displayGet(presentationDisplay, 'inverse gamma'));
    test.RGBimageGammaCorrected = test.RGBimageGammaCorrected / max(test.RGBimageGammaCorrected(:));

    %% Generate scene corresponding to the test stimulus on the presentation display
    theScene = sceneFromFile(test.RGBimageGammaCorrected,'rgb',background.luminance, presentationDisplay);
    theScene = sceneSet(theScene, 'h fov', fieldOfViewDegs);

    %% Verify that we have the desired cone contrast profiles
    % Get the emitted radiance image
    emittedRadianceImage = sceneGet(theScene, 'energy');
    
    %% Load the 2-deg Stockman cone fundamentals on a wavelength support matching the display
    coneFundamentals = ieReadSpectra(fullfile(isetbioDataPath,'human','stockman'), displayGet(presentationDisplay, 'wave'));
    
    %% Compute the LMS cone contrasts of the emitted radiance image
    test.achievedLMScontrastImage = computeLMScontrastImage(emittedRadianceImage, coneFundamentals);

    %% Visualize the L-cone, M-cone, S-cone components as well as the composite stimulus
    visualizeDisplayImage(1, test.RGBimage, test.LMSexcitationImage, presentationDisplay);
    
    %% Visualize the input cone contrasts and the achieved ones
    visualizeContrastImages(2, test.achievedLMScontrastImage, test.LMScontrastImage, fieldOfViewDegs);
end

function LMScontrastImage = computeLMScontrastImage(radianceImage, coneFundamentals)
    rowsNum = size(radianceImage,1);
    colsNum = size(radianceImage,2);
    wavelengthsNum = size(radianceImage, 3);
    radianceImage = reshape(radianceImage, [rowsNum*colsNum wavelengthsNum]);
        
    coneExcitationsImage = radianceImage * coneFundamentals;
    coneExcitationsBackground = coneExcitationsImage(1,:);

    LMScontrastImage = bsxfun(@times, ...
        bsxfun(@minus, coneExcitationsImage, coneExcitationsBackground), ...
        1./coneExcitationsBackground);
    
    LMScontrastImage = reshape(LMScontrastImage, [rowsNum colsNum 3]);
end

function LMScontrastImage = generateSpatialContrastPatterns(fieldOfViewDegs, spatialModulationParams)
    x = linspace(-fieldOfViewDegs/2,fieldOfViewDegs/2,256);
    [X,Y] = meshgrid(x);
    coneIDs = keys(spatialModulationParams);
    for coneIndex = 1:numel(coneIDs)
        p = spatialModulationParams(coneIDs{coneIndex});
        contrastPattern = cos(2*pi*(p.spatialFrequency(1)*(X-p.gaborPos(1)) + p.spatialFrequency(2)*(Y-p.gaborPos(2))));
        contrastEnvelope = exp(-((X-p.gaborPos(1))/p.gaborSigma(1)).^2) .* exp(-((Y-p.gaborPos(2))/p.gaborSigma(2)).^2);
        LMScontrastImage(:,:,coneIndex) = p.maxContrast * contrastPattern .* contrastEnvelope;
    end
end

function presentationDisplay = generatePresentationDisplay()
    % Generate presentation display
    presentationDisplay = displayCreate('LCD-Apple');
    bitDepth = 16;
    gammaT = 'nonlinear';
    if (strcmp(gammaT, 'linear'))
        N = 2^bitDepth;
        gTable = repmat(linspace(0, 1, N), 3, 1)';
    else
        x = linspace(0,1,2^bitDepth );
        gTable = x(:).^2.1;
        gTable = repmat(gTable, [1,3]);
    end
    presentationDisplay = displaySet(presentationDisplay, 'gTable', gTable);
end

function spatialModulationParams = getModulationParamsForStimType(stimType)
    spatialModulationParams = containers.Map();
    switch (stimType)
        case 'all-different'
            spatialModulationParams('L') = struct('maxContrast',  0.1, 'gaborPos', [0 0], 'gaborSigma', [Inf 0.6], 'spatialFrequency', [1.5 0]);
            spatialModulationParams('M') = struct('maxContrast', -0.1, 'gaborPos', [0 0], 'gaborSigma', [1.0 1.0], 'spatialFrequency', [0 1]);
            spatialModulationParams('S') = struct('maxContrast', -0.8, 'gaborPos', [0 0], 'gaborSigma', [1.2 1.2], 'spatialFrequency', [0.6 0]);
        case 'S-(L+M)'
            spatialModulationParams('L') = struct('maxContrast', -0.2, 'gaborPos', [0 0], 'gaborSigma', [1.0 1.0], 'spatialFrequency', [1 0]);
            spatialModulationParams('M') = struct('maxContrast', -0.2, 'gaborPos', [0 0], 'gaborSigma', [1.0 1.0], 'spatialFrequency', [1 0]);
            spatialModulationParams('S') = struct('maxContrast',  0.8, 'gaborPos', [0 0], 'gaborSigma', [1.0 1.0], 'spatialFrequency', [1 0]);
        case 'L-M'
            spatialModulationParams('L') = struct('maxContrast',  0.1, 'gaborPos', [0 0], 'gaborSigma', [1.0 1.0], 'spatialFrequency', [1 0]);
            spatialModulationParams('M') = struct('maxContrast', -0.1, 'gaborPos', [0 0], 'gaborSigma', [1.0 1.0], 'spatialFrequency', [1 0]);
            spatialModulationParams('S') = struct('maxContrast',  0.0, 'gaborPos', [0 0], 'gaborSigma', [1.0 1.0], 'spatialFrequency', [1 0]);
        case 'L+M'
            spatialModulationParams('L') = struct('maxContrast',  0.2, 'gaborPos', [0 0], 'gaborSigma', [1.0 1.0], 'spatialFrequency', [1 0]);
            spatialModulationParams('M') = struct('maxContrast',  0.2, 'gaborPos', [0 0], 'gaborSigma', [1.0 1.0], 'spatialFrequency', [1 0]);
            spatialModulationParams('S') = struct('maxContrast',  0.0, 'gaborPos', [0 0], 'gaborSigma', [1.0 1.0], 'spatialFrequency', [1 0]);
        case 'L+M+S'
            spatialModulationParams('L') = struct('maxContrast',  0.2, 'gaborPos', [0 0], 'gaborSigma', [1.0 1.0], 'spatialFrequency', [1 0]);
            spatialModulationParams('M') = struct('maxContrast',  0.2, 'gaborPos', [0 0], 'gaborSigma', [1.0 1.0], 'spatialFrequency', [1 0]);
            spatialModulationParams('S') = struct('maxContrast',  0.2, 'gaborPos', [0 0], 'gaborSigma', [1.0 1.0], 'spatialFrequency', [1 0]);
        otherwise
            error('Not know stimType: ''%s''.', stimType);
    end
end

function visualizeDisplayImage(figNo, RGBimage, LMSimage, presentationDisplay)

    hFig = figure(figNo);
    set(hFig, 'Position', [400 10 560 850]);
    % Generate axes in a [2x2] layout
    posVectors =  NicePlot.getSubPlotPosVectors(...
        'rowsNum', 2, 'colsNum', 2, ...
        'rightMargin', 0.00, ...
        'leftMargin', 0.01, ...
        'widthMargin', 0.01, ...
        'heightMargin', 0.05, ...
        'bottomMargin', 0.02, ...
        'topMargin', 0.03);
    
    % Get the bacground LMS excitations from the top-left pixel
    LMSbackground = LMSimage(1,1,:);
    
    % Extract the L-cone excitation component
    MSbackground = [0 LMSbackground(2) LMSbackground(3)];
    tmp = bsxfun(@plus, ...
                 bsxfun(@times, LMSimage, reshape([1 0 0], [1 1 3])), ...
                 reshape(MSbackground, [1 1 3]));
    LexcitationRGBimage = imageLinearTransform(tmp, inv(displayGet(presentationDisplay, 'rgb2lms')));
    
    % Extract the M-cone excitation component
    LSbackground = [LMSbackground(1) 0 LMSbackground(3)];
    tmp = bsxfun(@plus, ...
                 bsxfun(@times, LMSimage, reshape([0 1 0], [1 1 3])), ...
                 reshape(LSbackground, [1 1 3]));
    MexcitationRGBimage = imageLinearTransform(tmp, inv(displayGet(presentationDisplay, 'rgb2lms')));
    
    % Extract the S-cone excitation component
    LMbackground = [LMSbackground(1) LMSbackground(2) 0];
    tmp = bsxfun(@plus, ...
                 bsxfun(@times, LMSimage, reshape([0 0 1], [1 1 3])), ...
                 reshape(LMbackground, [1 1 3]));
    SexcitationRGBimage = imageLinearTransform(tmp, inv(displayGet(presentationDisplay, 'rgb2lms')));
    
    % Plot the L-cone component
    theCurrentAxes = subplot('position', posVectors(1,1).v);
    image(theCurrentAxes, lrgb2srgb(LexcitationRGBimage));
    axis(theCurrentAxes, 'square');
    set(theCurrentAxes, 'XTick', [], 'YTick', []);
    title(theCurrentAxes, 'L-cone stimulus component');
    
    % Plot the M-cone component
    theCurrentAxes = subplot('position', posVectors(1,2).v);
    image(theCurrentAxes, lrgb2srgb(MexcitationRGBimage));
    axis(theCurrentAxes, 'square');
    set(theCurrentAxes, 'XTick', [], 'YTick', []);
    title(theCurrentAxes, 'M-cone stimulus component');
    
    % Plot the S-cone component
    theCurrentAxes = subplot('position', posVectors(2,1).v);
    image(theCurrentAxes, lrgb2srgb(SexcitationRGBimage));
    axis(theCurrentAxes, 'square');
    set(theCurrentAxes, 'XTick', [], 'YTick', []);
    title(theCurrentAxes, 'S-cone stimulus component');
    
    % Plot the composite stimulus
    theCurrentAxes = subplot('position', posVectors(2,2).v);
    image(theCurrentAxes, lrgb2srgb(RGBimage));
    axis(theCurrentAxes, 'square')
    set(theCurrentAxes, 'XTick', [], 'YTick', []);
    title(theCurrentAxes, 'composite stimulus');
end

function visualizeContrastImages(figNo, outputLMScontrastImage, inputLMScontrastImage, fieldOfViewDegs)
    CLims = max([ ...
        max(abs(inputLMScontrastImage(:))) ...
        max(abs(outputLMScontrastImage(:))) ...
        ]) * [-1 1];
    
    hFig = figure(figNo); clf;
    set(hFig, 'Position', [1000 10 560 850]);
    % Generate axes in a [3x3] layout
    posVectors =  NicePlot.getSubPlotPosVectors(...
        'rowsNum', 3, 'colsNum', 3, ...
        'rightMargin', 0.00, ...
        'leftMargin', 0.01, ...
        'widthMargin', 0.01, ...
        'heightMargin', 0.05, ...
        'bottomMargin', 0.02, ...
        'topMargin', 0.03);
    
    theConeNames = {'L', 'M', 'S'};
    x = linspace(-fieldOfViewDegs/2,fieldOfViewDegs/2,size(inputLMScontrastImage,1));
    for k = 1:numel(theConeNames)
        theCurrentAxes = subplot('position', posVectors(1,k).v);
        imagesc(theCurrentAxes, x,x,squeeze(inputLMScontrastImage(:,:,k)), CLims); 
        axis(theCurrentAxes, 'square');
        set(theCurrentAxes, 'XTick', [], 'YTick', []);
        title(theCurrentAxes,sprintf('input %s-cone contrast', theConeNames{k}));
    end
    colormap(gray(1024))
    
    for k = 1:numel(theConeNames)
    	theCurrentAxes = subplot('position', posVectors(2,k).v);
        imagesc(theCurrentAxes, x, x, squeeze(outputLMScontrastImage(:,:,k)), CLims);
        axis(theCurrentAxes, 'square');
        set(theCurrentAxes, 'XTick', [], 'YTick', []);
        title(theCurrentAxes,sprintf('achieved %s-cone contrast', theConeNames{k}));
    end
    
    m = (size(outputLMScontrastImage,1))/2+1;
    for k = 1:numel(theConeNames)
        theCurrentAxes = subplot('position', posVectors(3,k).v);
        plot(theCurrentAxes,x, squeeze(inputLMScontrastImage(m,:,k)), 'k-', 'LineWidth', 3);
        hold(theCurrentAxes, 'on');
        plot(theCurrentAxes,x, squeeze(outputLMScontrastImage(m,:,k)), 'g--', 'LineWidth', 2);
        if (k == 3)
            legend(theCurrentAxes, {'input', 'achieved'}, 'Location', 'North', 'Orientation', 'horizontal');
        end
        axis(theCurrentAxes, 'square');
        if (k > 1)
            set(theCurrentAxes, 'YTickLabel', {});
        else
            ylabel(theCurrentAxes, 'contrast');
        end
        set(theCurrentAxes, 'XLim', [x(1) x(end)], 'XTick', [-5:1:5], 'YLim', CLims, 'YTick', -1:0.2:1, 'YLim', [-1 1]);
        box(theCurrentAxes, 'off');
        xlabel(theCurrentAxes, 'space (deg)');
    end
    
end
