function analyzeConeMosaic()

    subjectID = 4;
    mosaicEcc = [0 0];
    
    mosaicType = 'regular';
    analyze(mosaicType, mosaicEcc, subjectID, -4);
    
    mosaicType = 'ecc_varying';
    analyze(mosaicType, mosaicEcc, subjectID, -6);
    
end


function analyze(mosaicType, mosaicEcc, subjectID, minC)

    conesAcross = 75;
    cm = generateMosaic(mosaicType, conesAcross, mosaicEcc);

    % Generate psf
    
    wavefrontSpatialSamples = 901;
    psf = generatePSF(cm,subjectID,wavefrontSpatialSamples);
     
    % Generate PSF and coneMosaicImage
    fftSize = 2048;
    [psfImage, coneMosaicImage, pixelIndicesForConeAveraging, pixelConeIndexForConeAveraging, spatialSupportArcMin, ...
         Xhires,Yhires, NyquistCyclesPerDegreeFromConeSeparation] = generatePSFandConeImage(cm, psf, fftSize);
         
    % Compute cone mosaic sampling image
    coneMosaicSamplingImage = generateConeMosaicSamplingImage(coneMosaicImage, pixelIndicesForConeAveraging, ...
        pixelConeIndexForConeAveraging, fftSize);
    
    % Compute power spectrum of the cone mosaic
    [coneMosaicSpectrum, coneMosaicSpectrumSlice, spectralSupportCPD] = ...
             computeFourierPowerSpectrum(spatialSupportArcMin, coneMosaicSamplingImage, fftSize, NyquistCyclesPerDegreeFromConeSeparation);

    % Compute power spectrum of the PSF
    [psfSpectrum, psfSpectrumSlice, spectralSupportCPD] = ...
             computeFourierPowerSpectrum(spatialSupportArcMin, psfImage, fftSize, NyquistCyclesPerDegreeFromConeSeparation );
         
    
    % Examined spatial frequencies
    spatialFrequencyCPDexamined = 2:1:(ceil(2.5*NyquistCyclesPerDegreeFromConeSeparation));
    
    % Visualization setup
    rowsNum = 2;
    colsNum = 3;
    sv = NicePlot.getSubPlotPosVectors(...
           'colsNum', colsNum, ...
           'rowsNum', rowsNum, ...
           'heightMargin',  0.09, ...
           'widthMargin',    0.01, ...
           'leftMargin',     0.01, ...
           'rightMargin',    0.00, ...
           'bottomMargin',   0.04, ...
           'topMargin',      0.01); 
    
    videoFileName = sprintf('aliasingDemo_%sMosaic_ecc_%2.1degs%2.0f_subject',mosaicType, mosaicEcc(1), subjectID);
    videoOBJ = VideoWriter(videoFileName, 'MPEG-4');
    videoOBJ.FrameRate = 30;
    videoOBJ.Quality = 100;
    videoOBJ.open();
            
    hFigVideo = figure(1); clf;
    set(hFigVideo, 'Position', [10 10 1800 1360], 'Color', [1 1 1]);
    
    spatialRange = 13*[-1 1];
    sfRange = 145*[-1 1];
    
    % The mosaic with apertures shown as delta samples
    ax = subplot('Position', sv(2,1).v);
    imagesc(ax, spatialSupportArcMin, spatialSupportArcMin, coneMosaicSamplingImage);
    axis(ax, 'square'); axis(ax, 'xy');
    set(ax, 'XLim', spatialRange, 'YLim', spatialRange);
    set(ax, 'XTick', -60:5:60, 'YTick', []);
    set(ax, 'FontSize', 16);
    xlabel(ax, 'space (arc min)');
    colormap(ax, gray(1024))
    title(ax, sprintf('mosaic sampling grid (ecc = %2.1f, %2.1f degs)', mosaicEcc(1), mosaicEcc(2)));
    
    for iSF  = 1:numel(spatialFrequencyCPDexamined)
     
        spatialFrequencyCPD = spatialFrequencyCPDexamined(iSF);   
        
        % Make cone mosaic activation image and stimulus image
         [coneMosaicActivationImage, retinalStimulusImage] = ...
             generateConeMosaicActivationAndRetinalStimulusImage(coneMosaicImage, psfImage, spatialFrequencyCPD, Xhires/cm.micronsPerDegree, Yhires/cm.micronsPerDegree);
         
         % Make cone mosaic activation sampling image
         coneMosaicActivationSamplingImage = generateConeMosaicSamplingImage(coneMosaicActivationImage, pixelIndicesForConeAveraging, pixelConeIndexForConeAveraging, fftSize);
         
         % Normalize with respect to lowest frequency
         if (iSF == 1)
             maxActivationSamplingImage = max(abs(coneMosaicActivationSamplingImage(:)));
             maxActivationImage = max(abs(coneMosaicActivationImage(:)));
             maxRetinalStimulusImage = max(abs(retinalStimulusImage(:)));
         end
         coneMosaicActivationSamplingImage = coneMosaicActivationSamplingImage / maxActivationSamplingImage;
         coneMosaicActivationImage = coneMosaicActivationImage / maxActivationImage;
         retinalStimulusImage = retinalStimulusImage / maxRetinalStimulusImage;
         
         % Compute power spectrum of the stimulus
         [stimulusSpectrum, stimulusSpectrumSlice, spectralSupportCPD] = ...
             computeFourierPowerSpectrum(spatialSupportArcMin, retinalStimulusImage, fftSize, NyquistCyclesPerDegreeFromConeSeparation );
         
         % Compute power spectrum of the cone mosaic activation 
         [coneMosaicActivationSpectrum, coneMosaicActivationSpectrumSlice, spectralSupportCPD, reconstructedStimulusImage] = ...
             computeFourierPowerSpectrum(spatialSupportArcMin, coneMosaicActivationSamplingImage, fftSize, NyquistCyclesPerDegreeFromConeSeparation );
         
         
         % Normalize with respect to lowest frequency
         if (iSF == 1)
             maxStimulusSpectrum = max(stimulusSpectrum(:));
             maxActivationSpectrum = max(coneMosaicActivationSpectrum(:));
             maxReconstructedStimulusImage = max(abs(reconstructedStimulusImage(:)));
         end
         stimulusSpectrum = stimulusSpectrum / maxStimulusSpectrum;
         coneMosaicActivationSpectrum = coneMosaicActivationSpectrum / maxActivationSpectrum;
         reconstructedStimulusImage = reconstructedStimulusImage / maxReconstructedStimulusImage;

         ax = subplot('Position', sv(1,1).v);
         stimulusAssembly = reconstructedStimulusImage;
         midPoint = round(0.5*size(stimulusAssembly,1));
         stimulusAssembly(midPoint:end,:) = retinalStimulusImage(midPoint:end,:);
         imagesc(ax,spatialSupportArcMin, spatialSupportArcMin, stimulusAssembly, [-1 1]);
         hold(ax, 'on');
         plot(ax, [spatialSupportArcMin(1) spatialSupportArcMin(end)], [0 0], 'r--', 'LineWidth', 1.5);
         hold(ax, 'off');
         axis(ax, 'square'); axis(ax, 'xy');
         set(ax, 'XLim', spatialRange, 'YLim', spatialRange);
         set(ax, 'FontSize', 16);
         set(ax, 'XTick', -60:5:60, 'YTick', []);
         colormap(ax, gray(1024));
         title(ax, sprintf('input&reconstructed stimulus (%2.0f c/deg)',spatialFrequencyCPD));
         xlabel(ax, 'space (arc min)');
         
         % The cone mosaic cone activation (no spatial integration within apertures)
         ax = subplot('Position', sv(1,2).v);
         imagesc(ax, spatialSupportArcMin, spatialSupportArcMin, coneMosaicActivationImage, [-1 1]);
         set(ax, 'XLim', spatialRange, 'YLim', spatialRange);
         axis(ax, 'square'); axis(ax, 'xy');
         set(ax, 'XTick', -60:5:60, 'YTick', []);
         colormap(ax, gray(1024));
         set(ax, 'FontSize', 16);
         title(ax, 'retinal stimulus (within cone apertures)');
         xlabel(ax, 'space (arc min)');
         
         % The cone mosaic cone activation (spatial integration within apertures, and apertures as delta points)
         ax = subplot('Position', sv(2,2).v);
         imagesc(ax, spatialSupportArcMin, spatialSupportArcMin, coneMosaicActivationSamplingImage, [-1 1]);
         set(ax, 'XLim', spatialRange, 'YLim', spatialRange);
         axis(ax, 'square'); axis(ax, 'xy');
         set(ax, 'XTick', -60:5:60, 'YTick', []);
         set(ax, 'FontSize', 16);
         xlabel(ax, 'space (arc min)');
         colormap(ax, gray(1024))
         title(ax, 'mosaic sampling grid activation (integration within apertures)');
         
         % The spectrum of the stimulus
         ax = subplot('Position', sv(1,3).v);
         imagesc(ax, spectralSupportCPD, spectralSupportCPD, log10(stimulusSpectrum), [minC 0]);
         hold(ax, 'on');
         plot(ax,NyquistCyclesPerDegreeFromConeSeparation*cosd(0:5:360), ...
             NyquistCyclesPerDegreeFromConeSeparation*sind(0:5:360), 'c-', 'LineWidth', 1.5);
         hold(ax, 'off');
         axis(ax, 'square'); axis(ax, 'xy');
         set(ax, 'XLim', sfRange, 'Ylim', sfRange, 'XTick', -120:30:120, 'YTick', []);
         set(ax, 'FontSize', 16);
         colormap(ax, gray(1024))
         title(ax, 'stimulus spectrum');
         xlabel(ax,'spatial frequency (c/deg)');
         
         % The spectrum of cone mosaic activation (spatial integration within apertures, and apertures as delta points)
         ax = subplot('Position', sv(2,3).v);
         imagesc(ax,spectralSupportCPD, spectralSupportCPD, log10(coneMosaicActivationSpectrum), [minC 0]);
         hold(ax, 'on');
         plot(ax,NyquistCyclesPerDegreeFromConeSeparation*cosd(0:5:360), ...
             NyquistCyclesPerDegreeFromConeSeparation*sind(0:5:360), 'c-', 'LineWidth', 1.5);
         
         hold(ax, 'off');
         axis(ax, 'square'); axis(ax, 'xy');
         set(ax, 'XLim', sfRange, 'Ylim', sfRange, 'XTick', -120:30:120, 'YTick', []);
         set(ax, 'FontSize', 16);
         xlabel(ax,'spatial frequency (c/deg)');
         colormap(ax, gray(1024))
         title(ax, 'mosaic sampling grid activation spectrum');
         
         drawnow;
         videoOBJ.writeVideo(getframe(hFigVideo));
    end         
end


function [theImageSpectrum, theImageSpectrumSlice, spectralSupportCPD, reconstructionImage] = ...
    computeFourierPowerSpectrum(theSupportArcMin, theImage, fftSize, NyquistCyclesPerDegreeFromConeSeparation)

    % Zero padding
    n = fftSize-size(theImage,1);
    theImage = padarray(theImage,[n/2 n/2]);

    FT = fft2(theImage, fftSize, fftSize);
    theImageSpectrum = (abs(fftshift(FT))).^2;

    dxDegs = (theSupportArcMin(2)-theSupportArcMin(1))/60;
    maxSF = 1/(2*dxDegs);
    deltaSF = maxSF / (fftSize/2);
    spectralSupportCPD = ((-maxSF):deltaSF:maxSF-deltaSF);

    [X,Y] = meshgrid(spectralSupportCPD,spectralSupportCPD);
    R = sqrt(X.^2 + Y.^2);
    idx = find(R(:) < NyquistCyclesPerDegreeFromConeSeparation);
    reconstructionKernel = zeros(fftSize, fftSize);
    reconstructionKernel(idx) = 1;
    reconstructionKernel = ifftshift(reconstructionKernel);

    reconstructionFT = (abs(FT) .* reconstructionKernel) .* exp(1i * angle(FT));
    reconstructionImage = ifft2(reconstructionFT, fftSize, fftSize);
    %[max(real(reconstructionImage(:))) max(imag(reconstructionImage(:)))]

    reconstructionImage = real(reconstructionImage);
    %[min(reconstructionImage(:)) max(reconstructionImage(:))]

    % Generate 1D slice
    theImageSpectrumSlice = theImageSpectrum(fftSize/2,:);
end


function psfImageRotated = rotatedPSFimage(psfImage, fftSize)
    psfImageSpectrum = fftshift(abs(fft2(psfImage, fftSize,fftSize)));
    psfImageSpectrum = psfImageSpectrum / max(psfImageSpectrum(:));
    
    binaryImage = psfImageSpectrum;
    m1 = min(binaryImage(:));
    m2 = max(binaryImage(:));
    binaryImage = imbinarize((binaryImage - m1)/(m2-m1));
    s = regionprops(binaryImage,'orientation');
    psfRotation = s.Orientation;
    
    psfImageRotated = imrotate(psfImage, -psfRotation, 'bilinear', 'crop');
    psfImageRotated = psfImageRotated / max(psfImageRotated(:));
end

%
%
%
function cm = generateMosaic(type, conesAcross, mosaicEcc)
    % Mosic eccentricity and size
    whichEye = 'right eye';
   
    mosaicEccMicrons = 1e3 * RGCmodels.Watson.convert.rhoDegsToMMs(mosaicEcc);
    coneSpacingDegs = RGCmodels.Watson.compute.rfSpacingAtRetinalPositions(...
        'right eye', mosaicEccMicrons, 'cones', false);
    mosaicSizeDegs = coneSpacingDegs * conesAcross;
   
    % Generate mosaic centered at target eccentricity
    cm = cMosaic(...
        'sizeDegs', mosaicSizeDegs*[1 1], ...       % SIZE in degs
        'eccentricityDegs', mosaicEcc, ...          % ECC  in degs
        'whichEye', whichEye ...
        );

    if (strcmp(type, 'regular'))
        % Generate source mosaic, a @coneMosaicHex object
        meanConeSpacingMicrons = mean(cm.coneRFspacingsMicrons);
        coneApertureMicrons = 0.7 * meanConeSpacingMicrons;
        sourceMosaic = coneMosaicHex(7, ...
            'fovDegs', max(mosaicSizeDegs), ...
            'customInnerSegmentDiameter', coneApertureMicrons, ...
            'customLambda', meanConeSpacingMicrons, ...
            'micronsPerDegree', cm.micronsPerDegree);
        coneData = sourceMosaic.coneData();
        coneData.positions(:,1) = coneData.positions(:,1) + cm.eccentricityMicrons(1); 
        coneData.positions(:,2) = coneData.positions(:,2) + cm.eccentricityMicrons(2);
        % Generate equivalent @cMosaic object
        cm = cMosaic('coneData', coneData, ...
            'whichEye', whichEye);
    end
    
end

function psf = generatePSF(cm,subjectID,wavefrontSpatialSamples)

    [oiEnsemble, psfEnsemble] = ...
                cm.oiEnsembleGenerate(cm.eccentricityDegs, ...
                'zernikeDataBase', 'Polans2015', ...
                'subjectID', subjectID, ...
                'pupilDiameterMM', 3.0, ...
                'subtractCentralRefraction', false, ...
                'wavefrontSpatialSamples', wavefrontSpatialSamples, ...
                'zeroCenterPSF', true, ...
                'flipPSFUpsideDown', true);
    psf = psfEnsemble{1};       
end

function coneMosaicSamplingImage = generateConeMosaicSamplingImage(coneMosaicImage, pixelIndicesForConeAveraging, pixelConeIndexForConeAveraging, fftSize)
    
    coneMosaicSamplingImage = coneMosaicImage * 0;
    conesNum = numel(pixelConeIndexForConeAveraging);
    x = 1:fftSize;
    [Xhires, Yhires] = meshgrid(x,x);
    XY = [Xhires(:) Yhires(:)];
    
    for iCone = 1:conesNum
        idxCone = pixelConeIndexForConeAveraging(iCone);
        % average pixels within cone aperture
        pixelIndicesWithinConeAperture = pixelIndicesForConeAveraging{iCone};
        coneMosaicSamplingImage(XY(idxCone,2),XY(idxCone,1)) = mean(coneMosaicImage(pixelIndicesWithinConeAperture));
    end
    
end


function [psfImage, coneMosaicImage, pixelIndicesForConeAveraging, pixelConeIndexForConeAveraging, ...
          psfSupportArcMin, Xhires,Yhires, ...
          NyquistCyclesPerDegreeFromConeSeparation] = generatePSFandConeImage(cm, psf, fftSize)
      
     if isnan(min(psf.data(:))) || isnan(max(psf.data(:)))
         psfImage = [];
         coneMosaicImage = [];
         psfSupportArcMin = [];
         return;
     end
     
     % Make PSF image
     targetWavelength = 550;
     [~,idx] = min(abs(psf.supportWavelength-targetWavelength));
     psfImage = squeeze(psf.data(:,:,idx));
     
     % Rotate so that highest spatial frequency resolving axis is along the horizontal
     psfImage = rotatedPSFimage(psfImage, fftSize);
     
     % Transform psfSupport units from arc min to microns
     psfSupportMicrons = cm.micronsPerDegree * psf.supportX/60;
     [Xorig,Yorig] = meshgrid(psfSupportMicrons,psfSupportMicrons);
     
     % Interpolate to FFT samples
     interpolationSamples = fftSize;
     if (interpolationSamples < numel(psfSupportMicrons))
        error('Interpolation samples (%d) must be >= psf support length (%d)', ...
            interpolationSamples,numel(psfSupportMicrons));
     end
     
     psfSupportMicrons = linspace(psfSupportMicrons(1), psfSupportMicrons(end), interpolationSamples);
     [Xhires,Yhires] = meshgrid(psfSupportMicrons,psfSupportMicrons);
     psfImage = interp2(Xorig, Yorig, psfImage, Xhires, Yhires);
     psfImage = psfImage / max(psfImage(:));
     psfSupportArcMin = psfSupportMicrons / cm.micronsPerDegree * 60;
     
     % If there are no cones in this patch return an empty image
    if (isempty(cm.coneRFspacingsMicrons))
        coneMosaicImage = [];
        return;
    end
    
    % Centered cone positions
    conePositionsMicrons = bsxfun(@minus, cm.coneRFpositionsMicrons, cm.eccentricityMicrons);
    conesNum = size(conePositionsMicrons,1);
    % Cone apertures
    coneApertureMicrons = cMosaic.coneApertureToDiameterRatio * cm.coneRFspacingsMicrons;
    meanSeparationMicrons = median(cm.coneRFspacingsMicrons);
   
    meanSeparationDegs = meanSeparationMicrons / cm.micronsPerDegree;
    NyquistCyclesPerDegreeFromConeSeparation = 1/(2*meanSeparationDegs);
    
    % Generate cone mosaic image
    coneOutlineX = cosd(0:5:360);
    coneOutlineY = sind(0:5:360);
    coneMosaicImage = Xhires*0;

    XY = [Xhires(:) Yhires(:)];
    
    pixelIndicesForConeAveraging = cell(1, conesNum);
    pixelConeIndexForConeAveraging = zeros(1, conesNum);
    
    for iCone = 1:conesNum

        diffPos = bsxfun(@minus, XY,conePositionsMicrons(iCone,:));
        [~,idxCone] = min(sum(diffPos.^2,2));
        
        xx = conePositionsMicrons(iCone,1) + coneApertureMicrons(iCone)*0.5*coneOutlineX;
        yy = conePositionsMicrons(iCone,2) + coneApertureMicrons(iCone)*0.5*coneOutlineY;
        [in,on] = inpolygon(Xhires(:), Yhires(:), xx,yy);
        idx = find((in|on));
        
        % Add cone aperture
        tmp = Xhires*0;
        tmp(idx) = 1;
        coneMosaicImage = coneMosaicImage + tmp;

        % Keep track of averaging positions for this cone
        pixelIndicesForConeAveraging{iCone} = idx;
        pixelConeIndexForConeAveraging(iCone) = idxCone;
        
    end
end

function [coneMosaicActivationImage, stimulusImage] = ...
             generateConeMosaicActivationAndRetinalStimulusImage(coneMosaicImage, psfImage, stimulusSpatialFrequency, Xhires, Yhires)
      
    % Generate Gabor stimulus image
    orientationDegs = 0;
    stimulusImage = generateGaborImage(stimulusSpatialFrequency,orientationDegs, Xhires, Yhires);
    
    % Filter stimulus by PSF (this is a convolution in space, so multiplication in the frequency domain)
    fftSize = size(stimulusImage,1);
    stimulusPSFImage = ifftshift(ifft2(fft2(stimulusImage, fftSize,fftSize) .* fft2(psfImage, fftSize,fftSize), fftSize,fftSize));
    
    coneMosaicActivationImage = coneMosaicImage .* stimulusPSFImage;
end

function stimulusImage = generateGaborImage(spatialFrequencyCPD, orientationDegs, X, Y)
    sigma = max(X(:))/3;
    gaussianEnvelope = exp(-0.5*X.^2/sigma^2) .* exp(-0.5*Y.^2/sigma^2);
    stimulusImage = cos(2*pi*(X*spatialFrequencyCPD*cosd(orientationDegs) + Y*spatialFrequencyCPD*sind(orientationDegs)));
    stimulusImage = stimulusImage .* gaussianEnvelope;
end

