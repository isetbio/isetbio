function demoConeMosaicAliasing()

    % Get directory
    [directory,~] = fileparts(which(mfilename()));
    exportsDir = fullfile(directory, 'exports');
    
    
    whichEye = 'right eye';
     
    zernikeDataBase = 'Polans2015';  % choose between {'Polans2015', and 'Artal2012'}
    subjectID = 4;
    
    % Mosaic eccentricity
    mosaicEcc = [3 0];
    
    % How many cones across the mosaic to include
    conesAcross = 45;
    
    % Min spatial frequency to explore
    examinedSFmin = 4;
    
    % Spatial frequency support for visualized power spectrum
    sfRange = 165*[-1 1];
    
    % Spatial support for visualized point spread function
    spatialRange = 18*[-1 1];   
    
    % Perform analysis for a regular hex cone mosaic
    mosaicType = 'regular';
    analyze(whichEye, mosaicType, mosaicEcc, zernikeDataBase, subjectID, conesAcross , spatialRange, sfRange, examinedSFmin, -3.5, exportsDir);
    
    % Perform analysis for an eccentricity-varying cone mosaic
    mosaicType = 'ecc_varying';
    analyze(whichEye, mosaicType, mosaicEcc, zernikeDataBase, subjectID, conesAcross, spatialRange, sfRange, examinedSFmin, -5, exportsDir);
    
end


function analyze(whichEye, mosaicType, mosaicEcc, zernikeDataBase, subjectID, conesAcross, spatialRange, sfRange, examinedSFmin, minC, exportsDir)

    % Generate the cone mosaic
    cm = generateMosaic(whichEye, mosaicType, conesAcross, mosaicEcc);

    % Generate the psf at the mosaic's eccentricity
    wavefrontSpatialSamples = 901;
    psf = generatePSF(cm, zernikeDataBase, subjectID,wavefrontSpatialSamples);
     
    % Generate image of the PSF and of the coneMosaic
    fftSize = 2048;
    [psfImage, coneMosaicImage, pixelIndicesForConeAveraging, pixelConeIndexForConeAveraging, spatialSupportArcMin, ...
         Xhires,Yhires, NyquistCyclesPerDegreeFromConeSeparation] = generatePSFandConeImage(cm, psf, fftSize);
         
    % Generate image of the cone mosaic sampling with apertures depicted as delta functions
    coneMosaicSamplingImage = generateConeMosaicSamplingImage(coneMosaicImage, pixelIndicesForConeAveraging, ...
        pixelConeIndexForConeAveraging, fftSize);
    
    
    % Compute power spectrum of the cone mosaic
    [coneMosaicSpectrum, coneMosaicSpectrumSlice, spectralSupportCPD] = ...
             computeFourierPowerSpectrum(spatialSupportArcMin, coneMosaicSamplingImage, fftSize, NyquistCyclesPerDegreeFromConeSeparation);

    % Compute power spectrum of the PSF
    [psfSpectrum, psfSpectrumSlice, spectralSupportCPD] = ...
             computeFourierPowerSpectrum(spatialSupportArcMin, psfImage, fftSize, NyquistCyclesPerDegreeFromConeSeparation );
         
    
    % Examined spatial frequencies
    spatialFrequencyCPDexamined = examinedSFmin:1:(ceil(2.5*NyquistCyclesPerDegreeFromConeSeparation));
    
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
    
    % Video setup
    videoFileName = fullfile(exportsDir, sprintf('aliasingDemo_%sMosaic_ecc_%2.1degs_%2.1fdegs%2.0f_subject',mosaicType, mosaicEcc(1), mosaicEcc(2), subjectID));
    videoOBJ = VideoWriter(videoFileName, 'MPEG-4');
    videoOBJ.FrameRate = 30;
    videoOBJ.Quality = 100;
    videoOBJ.open();
            
    hFigVideo = figure(1); clf;
    set(hFigVideo, 'Position', [10 10 1800 1360], 'Color', [1 1 1]);
   
    
    % The mosaic with apertures shown as delta samples
    ax = subplot('Position', sv(2,1).v);
    imagesc(ax, spatialSupportArcMin, spatialSupportArcMin, psfImage); % coneMosaicSamplingImage);
    axis(ax, 'square'); axis(ax, 'xy');
    set(ax, 'XLim', spatialRange, 'YLim', spatialRange);
    set(ax, 'XTick', -60:5:60, 'YTick', []);
    set(ax, 'FontSize', 16);
    xlabel(ax, 'space (arc min)');
    colormap(ax, gray(1024))
    title(ax, sprintf('point spread function (ecc = %2.1f, %2.1f degs)', mosaicEcc(1), mosaicEcc(2)));
    
    maxActivationSamplingImage = 0;
    maxActivationImage = 0;
    maxRetinalStimulusEnergy = 0;
    maxStimulusSpectrum = 0;
    maxActivationSpectrum = 0;
    maxReconstructedStimulusEnergy = 0;
             
    for iSF  = 1:numel(spatialFrequencyCPDexamined)
     
        spatialFrequencyCPD = spatialFrequencyCPDexamined(iSF);   
        
        % Make cone mosaic activation image and stimulus image
         [coneMosaicActivationImage, retinalStimulusImage] = ...
             generateConeMosaicActivationAndRetinalStimulusImage(coneMosaicImage, psfImage, spatialFrequencyCPD, Xhires/cm.micronsPerDegree, Yhires/cm.micronsPerDegree);
         
         % Make cone mosaic activation sampling image
         coneMosaicActivationSamplingImage = generateConeMosaicSamplingImage(coneMosaicActivationImage, pixelIndicesForConeAveraging, pixelConeIndexForConeAveraging, fftSize);
         
         % Normalize
         maxActivationSamplingImageTmp = max(abs(coneMosaicActivationSamplingImage(:)));
         maxActivationImageTmp = max(abs(coneMosaicActivationImage(:)));
         maxRetinalStimulusEnergyTmp = sqrt(sum((retinalStimulusImage(:).^2)));
         maxActivationSamplingImage = max([maxActivationSamplingImageTmp maxActivationSamplingImage]);
         maxActivationImage = max([maxActivationImageTmp maxActivationImage]);
         maxRetinalStimulusEnergy = max([maxRetinalStimulusEnergyTmp maxRetinalStimulusEnergy]);

         coneMosaicActivationSamplingImage = coneMosaicActivationSamplingImage / maxActivationSamplingImage;
         coneMosaicActivationImage = coneMosaicActivationImage / maxActivationImage;
         
         % Scaling factor for retinal stimulus
         scalingFactor = sqrt(sum((retinalStimulusImage(:).^2)))/maxRetinalStimulusEnergy;
         retinalStimulusImage = retinalStimulusImage/max(abs(retinalStimulusImage(:))) * scalingFactor;
         
         % Compute power spectrum of the stimulus
         [stimulusSpectrum, stimulusSpectrumSlice, spectralSupportCPD] = ...
             computeFourierPowerSpectrum(spatialSupportArcMin, retinalStimulusImage, fftSize, NyquistCyclesPerDegreeFromConeSeparation );
         
         % Compute power spectrum of the cone mosaic activation 
         [coneMosaicActivationSpectrum, coneMosaicActivationSpectrumSlice, spectralSupportCPD, reconstructedStimulusImage] = ...
             computeFourierPowerSpectrum(spatialSupportArcMin, coneMosaicActivationSamplingImage, fftSize, NyquistCyclesPerDegreeFromConeSeparation );
         
         
         % Normalize with respect to lowest frequency
         maxStimulusSpectrumTmp = max(stimulusSpectrum(:));
         maxActivationSpectrumTmp = max(coneMosaicActivationSpectrum(:));
         maxReconstructedStimulusEnergyTmp = sqrt(sum((reconstructedStimulusImage(:).^2)));

         maxStimulusSpectrum = max([maxStimulusSpectrumTmp maxStimulusSpectrum]);
         maxActivationSpectrum = max([maxActivationSpectrumTmp maxActivationSpectrum]);
         maxReconstructedStimulusEnergy = max([maxReconstructedStimulusEnergyTmp maxReconstructedStimulusEnergy]);
         stimulusSpectrum = stimulusSpectrum / maxStimulusSpectrum;
         coneMosaicActivationSpectrum = coneMosaicActivationSpectrum / maxActivationSpectrum;

         % Scaling factor for reconstructed retinal stimulus
         scalingFactor = sqrt(sum((reconstructedStimulusImage(:).^2)))/maxReconstructedStimulusEnergy;
         reconstructedStimulusImage = reconstructedStimulusImage/max(abs(reconstructedStimulusImage(:))) * scalingFactor;
         
         ax = subplot('Position', sv(1,1).v);
         stimulusAssembly = reconstructedStimulusImage;
         midPoint = round(0.5*size(stimulusAssembly,1));
         stimulusAssembly(midPoint+1:end,:) = retinalStimulusImage(midPoint+1:end,:);
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

    % Compute power spectrum
    FT = fft2(theImage, fftSize, fftSize);
    theImageSpectrum = (abs(fftshift(FT))).^2;

    % Compute power spectrum spatial frequency support
    dxDegs = (theSupportArcMin(2)-theSupportArcMin(1))/60;
    maxSF = 1/(2*dxDegs);
    deltaSF = maxSF / (fftSize/2);
    spectralSupportCPD = ((-maxSF):deltaSF:maxSF-deltaSF);

    % Compute reconstruction kernel (cylinder within the Nyquist frequency)
    [X,Y] = meshgrid(spectralSupportCPD,spectralSupportCPD);
    R = sqrt(X.^2 + Y.^2);
    idx = find(R(:) < NyquistCyclesPerDegreeFromConeSeparation);
    reconstructionKernel = zeros(fftSize, fftSize);
    reconstructionKernel(idx) = 1;
   
    % Compute FT of reconstruction kernel
    reconstructionKernelFT = ifftshift(reconstructionKernel);
    
    % Perform filtering in frequency domain
    reconstructionFT = (abs(FT) .* reconstructionKernelFT) .* exp(1i * angle(FT));
    
    % Back to space domain
    reconstructionImage = ifft2(reconstructionFT, fftSize, fftSize);
    
    % Keep real part
    reconstructionImage = real(reconstructionImage);

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
function cm = generateMosaic(whichEye, type, conesAcross, mosaicEcc)
    % Mosic eccentricity and size
   
    mosaicEccMicrons = 1e3 * RGCmodels.Watson.convert.rhoDegsToMMs(mosaicEcc);
    coneSpacingDegs = RGCmodels.Watson.compute.rfSpacingAtRetinalPositions(...
        whichEye, mosaicEccMicrons, 'cones', false);
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

function psf = generatePSF(cm, zernikeDataBase, subjectID,wavefrontSpatialSamples)

    [oiEnsemble, psfEnsemble] = ...
                cm.oiEnsembleGenerate(cm.eccentricityDegs, ...
                'zernikeDataBase', zernikeDataBase, ...
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
        if (~isempty(pixelIndicesWithinConeAperture))
            coneMosaicSamplingImage(XY(idxCone,2),XY(idxCone,1)) = mean(coneMosaicImage(pixelIndicesWithinConeAperture));
        end
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

function [coneMosaicActivationImage, stimulusPSFImage] = ...
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

