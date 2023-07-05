function testOTF
    s_initISET;

    h1 = figure(1);
    set(h1, 'Position', [100 100 650 760]);
    clf;
    
    h2 = figure(2);
    set(h2, 'Position', [200 200 650 760]);
    clf;
    
    
    % Pupil diameters to test
    pupilDiametersInMillimeters = (2:0.5:6.5);
    %pupilDiametersInMillimeters = (4.0:0.5:6);
    
    for pupilSizeIndex = 1:numel(pupilDiametersInMillimeters)
        
        %% Retrieve examined pupil radius 
        pupilRadiusInMeters = pupilDiametersInMillimeters(pupilSizeIndex)/2.0/1000.0;
        
        %% Create human optics with given pupil radius
        optics = opticsCreate('human', pupilRadiusInMeters);
        
        %% Initialize optical image with above optics
        oi = oiCreate('human');
        oi = oiSet(oi, 'optics', optics);
        
        %% Compute optical image for given scene
        scene = sceneCreate('line d65');
        %% Make the scene angular size = 1 deg and place it at a distance = 1.0 m
        sceneAngularSizeInDeg = 2.0;
        sceneDistanceInMeters = 1.0;
        scene = sceneSet(scene,'wangular', sceneAngularSizeInDeg);
        scene = sceneSet(scene,'distance', sceneDistanceInMeters);
        
        %% Compute optical image
        oi = oiCompute(scene,oi); 
        
        %% Compute RGB rendition of optical image
        opticalRGBImage = oiGet(oi, 'rgb image');
        sceneRGBImage   = sceneGet(scene, 'rgb image');
        compositeRGBimage = generateCompositeImage(sceneRGBImage, opticalRGBImage);
        
        %% Retrieve the full OTF
        optics = oiGet(oi, 'optics');
        OTF    = abs(opticsGet(optics,'otf data'));
        
        %% Retrieve the wavelength axis
        OTFwavelengths = opticsGet(optics,'otf wave');
        
        %% Retrieve the spatial frequency support. This is in cycles/micron
        OTFsupport = opticsGet(optics,'otf support', 'um');
        
        otf_sfXInCyclesPerMicron = OTFsupport{1};
        otf_sfYInCyclesPerMicron = OTFsupport{2};
        
        %% Convert to cycles/deg.
        % In human retina, 1 deg of visual angle is about 288 microns
        micronsPerDegee = 288;
        otf_sfX = otf_sfXInCyclesPerMicron * micronsPerDegee;
        otf_sfY = otf_sfYInCyclesPerMicron * micronsPerDegee;
        
        %% Get the 2D slice at 550 nm
        [~,waveIndex]   = min(abs(OTFwavelengths - 550));
        
        %% Shift (0,0) to origin
        OTF550 = fftshift(squeeze(OTF(:,:,waveIndex)));
        
        %% Get a 2D slice through origin
        [~, sfIndex] = min(abs(otf_sfY - 0));
        OTFslice = squeeze(OTF550(sfIndex,:));
        
        %% Generate plots
        plotWidth = 0.3;
        plotHeight = 0.85/(numel(pupilDiametersInMillimeters)/2);
        margin     = 0.021;
        if (pupilSizeIndex <= numel(pupilDiametersInMillimeters)/2)
            figure(h1);
            subplotRow = pupilSizeIndex;
        else
            figure(h2);
            subplotRow = pupilSizeIndex-numel(pupilDiametersInMillimeters)/2;
        end
        
        %% Plot the 2D OTF
        subplot('Position', [0.03 1+margin/2-subplotRow*(plotHeight+margin) plotWidth plotHeight]);
        imagesc(otf_sfX, otf_sfY, OTF550);
        axis 'image'
        axis 'xy'
        if (pupilSizeIndex == numel(pupilDiametersInMillimeters))
           xlabel('cycles/deg');
        else
           xlabel(''); 
        end
        set(gca, 'XLim', [-60 60], 'YLim', [-60 60], 'XTick', [-60:60:60], 'YTick', [-60:60:60]);
        colormap(gray(256));
        
        %% Plot the 1D OTF slice
        subplot('Position', [0.06+plotWidth 1+margin/2-subplotRow*(plotHeight+margin) plotWidth plotHeight]);
        indices = find(otf_sfX >= 0);
        OTFsfX = otf_sfX(indices)+0.1;
        plot(OTFsfX, OTFslice(indices), 'rs-');
        if (pupilSizeIndex == numel(pupilDiametersInMillimeters))
           xlabel('cycles/deg');
        else
           xlabel(''); 
        end
        text(0.12, 0.005, sprintf('Pupil:%2.1f mm', pupilDiametersInMillimeters(pupilSizeIndex)), 'FontSize', 12, 'FontWeight', 'bold');
        set(gca, 'XLim', [0.1 100], 'YLim', [0.001 1]);
        set(gca, 'XScale', 'log', 'YScale', 'log', 'XTick', [0.1 1 2 5 10 20 50 100], 'YTick', [0.001 0.002 0.005 0.01 0.02 0.05 0.10 0.20 0.50 1.0]);
        set(gca, 'XTickLabel', [0.1 1 2 5 10 20 50 100], 'YTickLabel',  [0.001 0.002 0.005 0.01 0.02 0.05 0.10 0.20 0.50 1.0]);
        box on;
        grid on;
        colormap(gray(256));
        
        %% Plot the RGB rendered optical image
        subplot('Position', [0.09+2*plotWidth 1+margin/2-subplotRow*(plotHeight+margin) plotWidth plotHeight]);
        halfRange = size(compositeRGBimage,1)/2-10;
        rowRange = size(compositeRGBimage,1)/2 + [-halfRange:halfRange];
        colRange = size(compositeRGBimage,2)/2 + [-halfRange:halfRange];
        imshow(compositeRGBimage(rowRange, colRange, :));
        axis 'image'
        
        drawnow; 
    end
end

function compositeRGBimage = generateCompositeImage(sceneRGBImage, opticalRGBImage)
    compositeRGBimage = opticalRGBImage;
    dstRows = size(opticalRGBImage,1);
    dstCols = size(opticalRGBImage,2);
    srcRows = size(sceneRGBImage, 1);
    srcCols = size(sceneRGBImage, 2);
    rowMargin = (dstRows-srcRows)/2;
    colMargin = (dstCols-srcCols)/2;
    compositeRGBimage(1:dstRows/2,:,:) = 0;
    compositeRGBimage(rowMargin+(1:srcRows/2), colMargin+(1:srcCols/2),:) = sceneRGBImage(1:srcRows/2,1:srcCols/2,:);
end
