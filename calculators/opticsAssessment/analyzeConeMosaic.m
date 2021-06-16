function analyzeConeMosaic()

    % Mosic eccentricity and size
    mosaicEcc = [-7 -25];
    mosaicSizeDegs = 0.5*max([0.5 0.75*sqrt(sqrt(sum(mosaicEcc.^2,2)))])

    whichEye = 'right eye';
    
    % Generate mosaic centered at target eccentricity
    cm = cMosaic(...
        'sizeDegs', mosaicSizeDegs*[1 1], ...       % SIZE in degs
        'eccentricityDegs', mosaicEcc, ...          % ECC in degs
        'whichEye', whichEye ...
        );
    
    subjectID = 6;
    fftSize = 2048;
        
    
    [oiEnsemble, psfEnsemble] = ...
                cm.oiEnsembleGenerate(mosaicEcc, ...
                'zernikeDataBase', 'Polans2015', ...
                'subjectID', subjectID, ...
                'pupilDiameterMM', 3.0, ...
                'subtractCentralRefraction', false, ...
                'wavefrontSpatialSamples', 1001, ...
                'zeroCenterPSF', true, ...
                'flipPSFUpsideDown', true);
    psf = psfEnsemble{1};       
    
    
         
   % Make PSF and cone mosaic image
   coneApertureFraction = 0.01;
   
   includeOpticalFiltering = true;
   includeApertureFiltering =~true;
   if (coneApertureFraction <= 0.1)
       includeApertureFiltering =~true;
   end
   
   [psfImage, coneMosaicImage, psfSupportArcMin, Xhires,Yhires, NyquistCyclesPerDegreeFromConeSeparation] = ...
             makePSFandConeImage(cm, psf, coneApertureFraction, fftSize);
         
   if (~includeOpticalFiltering)
      psfImage = psfImage*0+1; 
   end
   
   YellotCircle.x = 2*NyquistCyclesPerDegreeFromConeSeparation * cosd(0:1:360);
   YellotCircle.y = 2*NyquistCyclesPerDegreeFromConeSeparation * sind(0:1:360);
         
         
    rowsNum = 3;
    colsNum = 3;
    sv = NicePlot.getSubPlotPosVectors(...
           'colsNum', colsNum, ...
           'rowsNum', rowsNum, ...
           'heightMargin',  0.07, ...
           'widthMargin',    0.01, ...
           'leftMargin',     0.02, ...
           'rightMargin',    0.00, ...
           'bottomMargin',   0.04, ...
           'topMargin',      0.03); 
       

     spatialFrequencyCPDexamined = 1:1:60;
     
     
     if (includeOpticalFiltering)
         if (includeApertureFiltering)
            videoOBJ = VideoWriter(sprintf('AliasingSubject_%2.0f_Ecc_%2.0f_%2.0fdegs.mp4', subjectID, mosaicEcc(1), mosaicEcc(2)), 'MPEG-4');
         else
             videoOBJ = VideoWriter(sprintf('AliasingSubjectNoAperture_%2.0f_Ecc_%2.0f_%2.0fdegs.mp4', subjectID, mosaicEcc(1), mosaicEcc(2)), 'MPEG-4');
         end
     else
         if (includeApertureFiltering)
            videoOBJ = VideoWriter(sprintf('AliasingNoOptics_Ecc_%2.0f_%2.0fdegs.mp4', mosaicEcc(1). mosaicEcc(2)), 'MPEG-4');
         else
             videoOBJ = VideoWriter(sprintf('AliasingNoOpticsNoAperture_Ecc_%2.0f_%2.0fdegs.mp4', mosaicEcc(1), mosaicEcc(2)), 'MPEG-4');
         end
     end
     
     videoOBJ.FrameRate = 30;
     videoOBJ.Quality = 100;
     videoOBJ.open();
     
     powerSpectrumRange = [-4 0];
     
     maxActivationSpectrum = [];
     indicesList = [];
     
     
     hFig = figure(1); clf;
     set(hFig, 'Position', [10 10 1500 1500]);
         
     r = sqrt(sum(mosaicEcc.^2,2));
     if (r > 1)
        maxVisualizedSF = 180/sqrt(r);
     end
     
     if (maxVisualizedSF >= 120)
         SFtickInterval = 60;
     elseif (maxVisualizedSF >= 90)
          SFtickInterval = 30;
     else
         SFtickInterval = 10;
     end
     
     for iSF  = 1:numel(spatialFrequencyCPDexamined)
     
         spatialFrequencyCPD = spatialFrequencyCPDexamined(iSF);
         
         % Make cone mosaic activation images
         [coneMosaicActivationImage, stimulusImage, indicesList] = ...
             makeConeMosaicActivationImage(cm, coneApertureFraction, psfImage, spatialFrequencyCPD,  Xhires,Yhires, indicesList);


         % Compute power spectrum of the stimulus
         [stimulusSpectrum, stimulusSpectrumSlice, spectralSupportCPD] = ...
             computeFourierPowerSpectrum(psfSupportArcMin, stimulusImage, [], fftSize );
         
         % Compute power spectrum of the cone mosaic
         [coneMosaicSpectrum, coneMosaicSpectrumSlice, spectralSupportCPD] = ...
             computeFourierPowerSpectrum(psfSupportArcMin, coneMosaicImage, [], fftSize );

         % Compute power spectrum of the PSF
         [psfSpectrum, psfSpectrumSlice, spectralSupportCPD] = ...
             computeFourierPowerSpectrum(psfSupportArcMin, psfImage, [], fftSize );
         
        
         % Compute power spectrum of the cone mosaic activation by the optical image of the sinewave stimulus
         [coneMosaicActivationSpectrum, coneMosaicActivationSpectrumSlice, ~, maxActivationSpectrum] = ...
                computeFourierPowerSpectrum(psfSupportArcMin, coneMosaicActivationImage, maxActivationSpectrum, fftSize);


         
         if (iSF == 1)
             ax = subplot('Position', sv(1,1).v);
             cm.visualize('figureHandle', hFig, 'axesHandle', ax, ...
                 'plotTitle', sprintf('cone mosaic \necc x=%2.1f, y=%2.1f', mosaicEcc(1), mosaicEcc(2)));
             hold(ax, 'on');
             semiTransparentContourPlot(ax, psfSupportArcMin/60+mosaicEcc(1), psfSupportArcMin/60+mosaicEcc(2), psfImage, 0.1:0.1:0.9, ...
                brewermap(1024, 'blues'), 0.3, [0 0 0]);
             hold(ax, 'off');
             
             ax = subplot('Position', sv(2,1).v);
             imagesc(ax,spectralSupportCPD, spectralSupportCPD, log10(coneMosaicSpectrum));
             hold(ax, 'on');
             % The Nyquist circle
             plot(ax, YellotCircle.x/2, YellotCircle.y/2, 'k-', 'LineWidth', 3.0);
             plot(ax, YellotCircle.x/2, YellotCircle.y/2, 'g-', 'LineWidth', 1.5);
             % The Yellot circle
             plot(ax, YellotCircle.x, YellotCircle.y, 'k-', 'LineWidth', 3.0);
             plot(ax, YellotCircle.x, YellotCircle.y, 'y-', 'LineWidth', 1.5);
             hold(ax, 'off');
             axis(ax, 'image'); axis 'xy';
             set(ax, 'XLim', maxVisualizedSF*[-1 1], 'YLim', maxVisualizedSF*[-1 1], ...
                 'XTick', [-240:SFtickInterval:240], 'YTick', [-240:SFtickInterval:240], ...
                 'CLim', powerSpectrumRange, 'FontSize', 16);
             grid on;
             xlabel(ax,'spatial frequency (c/deg)'); 
             ylabel('spatial frequency (c/deg)'); 
             colormap(ax,gray(1024));
             title(ax, 'cone mosaic spectrum');
             
             
             ax = subplot('Position', sv(3,1).v);
             imagesc(ax,spectralSupportCPD, spectralSupportCPD, log10(psfSpectrum));
             hold(ax, 'on');
             % The Nyquist circle
             plot(ax, YellotCircle.x/2, YellotCircle.y/2, 'k-', 'LineWidth', 3.0);
             plot(ax, YellotCircle.x/2, YellotCircle.y/2, 'g-', 'LineWidth', 1.5);
             % The Yellot circle
             plot(ax, YellotCircle.x, YellotCircle.y, 'k-', 'LineWidth', 3.0);
             plot(ax, YellotCircle.x, YellotCircle.y, 'y-', 'LineWidth', 1.5);
             hold(ax, 'off');
             axis(ax, 'image'); axis 'xy';
             set(ax, 'XLim', maxVisualizedSF*[-1 1], 'YLim', maxVisualizedSF*[-1 1], ...
                 'XTick', [-240:SFtickInterval:240], 'YTick', [-240:SFtickInterval:240], ...
                 'CLim', powerSpectrumRange, 'FontSize', 16);
             grid on;
             xlabel('spatial frequency (c/deg)'); 
             ylabel('spatial frequency (c/deg)'); 
             colormap(ax,gray(1024));
             title(ax, 'PSF spectrum');
             
             drawnow;
         end
         
         ax = subplot('Position', sv(1,2).v);
         imagesc(ax,psfSupportArcMin/60, psfSupportArcMin/60, stimulusImage);
         xlabel(ax, 'space (degrees)');
         axis(ax,'square'); axis 'xy';
         set(ax, 'XLim', mosaicSizeDegs*0.5*[-1 1], 'YLim', mosaicSizeDegs*0.5*[-1 1], 'FontSize', 16);
         if (mosaicSizeDegs < 0.3)
            set(ax, 'XTick', -0.2:0.1:0.2, 'YTick', -0.2:0.1:0.2);
         elseif (mosaicSizeDegs < 1)
             set(ax, 'XTick', -0.5:0.25:0.5, 'YTick', -0.5:0.25:0.5);
         else
             set(ax, 'XTick', -1:0.5:1, 'YTick', -1:0.5:1);
         end
         
         colormap(ax,gray(1024));
         title(ax, sprintf('stimulus frequency: %2.0f c/deg', spatialFrequencyCPD));
         
         
         
         ax = subplot('Position', sv(2,2).v);
         imagesc(ax,spectralSupportCPD, spectralSupportCPD, log10(coneMosaicActivationSpectrum));
         hold(ax, 'on');
         % The Nyquist circle
         plot(ax, YellotCircle.x/2, YellotCircle.y/2, 'k-', 'LineWidth', 3.0);
         plot(ax, YellotCircle.x/2, YellotCircle.y/2, 'g-', 'LineWidth', 1.5);
         % The Yellot circle
         plot(ax, YellotCircle.x, YellotCircle.y, 'k-', 'LineWidth', 3.0);
         plot(ax, YellotCircle.x, YellotCircle.y, 'y-', 'LineWidth', 1.5);
         hold(ax, 'off');
         axis(ax, 'image'); axis 'xy';
         set(ax, 'XLim', maxVisualizedSF*[-1 1], 'YLim', maxVisualizedSF*[-1 1], ...
             'XTick', [-240:SFtickInterval:240], 'YTick', [-240:SFtickInterval:240], ...
             'CLim', powerSpectrumRange , 'FontSize', 16);
         grid on;
         
         colormap(ax,gray(1024));
         title(ax, 'cone mosaic activation spectrum');
         
 
         ax = subplot('Position', sv(2,3).v);
         imagesc(ax,spectralSupportCPD, spectralSupportCPD, log10(stimulusSpectrum));
         hold(ax, 'on');
         % The Nyquist circle
         plot(ax, YellotCircle.x/2, YellotCircle.y/2, 'k-', 'LineWidth', 3.0);
         plot(ax, YellotCircle.x/2, YellotCircle.y/2, 'g-', 'LineWidth', 1.5);
         % The Yellot circle
         plot(ax, YellotCircle.x, YellotCircle.y, 'k-', 'LineWidth', 3.0);
         plot(ax, YellotCircle.x, YellotCircle.y, 'y-', 'LineWidth', 1.5);
         hold(ax, 'off');
         axis(ax, 'image'); axis 'xy';
         set(ax, 'XLim', maxVisualizedSF*[-1 1], 'YLim', maxVisualizedSF*[-1 1], ...
             'XTick', [-240:SFtickInterval:240], 'YTick', [-240:SFtickInterval:240], ...
             'CLim', powerSpectrumRange , 'FontSize', 16);
         grid on;
         xlabel(ax,'spatial frequency (c/deg)'); 
         colormap(ax,gray(1024));
         title(ax, 'stimulus spectrum');
         
         
         ax = subplot('Position', sv(3,2).v);
         % They Nyquist
         plot(ax, NyquistCyclesPerDegreeFromConeSeparation*[1 1], 10.^powerSpectrumRange, 'k-', 'LineWidth', 3);
         hold(ax, 'on');
         plot(ax,-NyquistCyclesPerDegreeFromConeSeparation*[1 1], 10.^powerSpectrumRange, 'k-', 'LineWidth', 3);
         plot(ax, NyquistCyclesPerDegreeFromConeSeparation*[1 1], 10.^powerSpectrumRange, 'g-', 'LineWidth', 1.5);
         plot(ax,-NyquistCyclesPerDegreeFromConeSeparation*[1 1], 10.^powerSpectrumRange, 'g-', 'LineWidth', 1.5);
         % The Yellot circle
         plot(ax, 2*NyquistCyclesPerDegreeFromConeSeparation*[1 1], 10.^powerSpectrumRange, 'k-', 'LineWidth', 3);
         plot(ax,-2*NyquistCyclesPerDegreeFromConeSeparation*[1 1], 10.^powerSpectrumRange, 'k-', 'LineWidth', 3);
         plot(ax, 2*NyquistCyclesPerDegreeFromConeSeparation*[1 1], 10.^powerSpectrumRange, 'y-', 'LineWidth', 1.5);
         plot(ax,-2*NyquistCyclesPerDegreeFromConeSeparation*[1 1], 10.^powerSpectrumRange, 'y-', 'LineWidth', 1.5);
         
         area(ax, spectralSupportCPD, psfSpectrumSlice, 10.^powerSpectrumRange(1), ...
             'FaceColor', [0.1 0.7 0.8], 'FaceAlpha', 0.5, 'EdgeColor', [0 .7 0.8], 'LineWidth', 0.5);
         area(ax, spectralSupportCPD, coneMosaicActivationSpectrumSlice, 10.^powerSpectrumRange(1), ...
             'FaceColor', [1.0 0.1 0.3], 'FaceAlpha', 0.5, 'EdgeColor', [1 0 0], 'LineWidth', 0.5);
         hold(ax, 'off');
         axis 'square'
         set(ax, 'XLim', maxVisualizedSF*[-1 1], 'XTick', -240:60:240, 'YLim', 10.^powerSpectrumRange);
         set(ax, ...
             'YTick', [0.00001  0.00003      0.0001  0.0003  0.001  0.003 0.01   0.03 0.1   0.3 1], ...
             'YTickLabel', {'1e-5', '3e-5', '1e-4', '3e-4', '1e-3', '3e-3', '1e-2',  '3e-2', '1e-1',  '3e-1', '1'} );
         set(ax, 'XScale', 'linear', 'YScale', 'log', 'FontSize', 16);
         set(ax, 'XLim', maxVisualizedSF*[-1 1], 'XTick', -240:SFtickInterval:240, ...
             'YLim', 10.^powerSpectrumRange);
         grid(ax,'on');
         ylabel(ax, 'power');
         xlabel(ax,'spatial frequency (c/deg)'); 
         
         drawnow;
         videoOBJ.writeVideo(getframe(hFig));
     end
     
     videoOBJ.close();
     
end

function [theImageSpectrum, theImageSpectrumSlice, spectralSupportCPD, maxActivationSpectrum] = ...
    computeFourierPowerSpectrum(theSupportArcMin, theImage, maxActivationSpectrum, fftSize)

    % Zero padding

     n = fftSize-size(theImage,1);
     theImage = padarray(theImage,[n/2 n/2]);
    
    % Window it
%     fftSize = size(theImage,1);
%     x = 1:fftSize;
%     x = x - mean(x);
%     [X,Y] = meshgrid(x,x);
%     
%     w = hamming(fftSize);
%     w = w*w';
%      

    FT = fftshift(fft2(theImage, fftSize, fftSize));
    F = abs(FT);
    F = F .^ 2;
    
    if (isempty(maxActivationSpectrum))
        maxActivationSpectrum = max(F(:));
    end
    
    theImageSpectrum = F / maxActivationSpectrum;
    
    
    dxDegs = (theSupportArcMin(2)-theSupportArcMin(1))/60;
    samplingFrequency = 1/dxDegs;
    nyquistFrequency = samplingFrequency/2;
    dF = samplingFrequency/fftSize;
    spectralSupportCPD = (-nyquistFrequency):dF:(nyquistFrequency-dF);
    
   
    % Generate 1D slice
    
     theImageSpectrumSlice = theImageSpectrum(fftSize/2,:);
    
%     deltaAngle = 5;
%     angles = 0:deltaAngle:(360-deltaAngle);
%     for iAngle = 1:numel(angles)
%         tmp = imrotate(theImageSpectrum, iAngle, 'bilinear', 'crop');
%         slice = slice + tmp(fftSize/2,:);
%     end
% 
%     theImageSpectrumSlice = slice / numel(angles);
    
    
end



function sinwaveImage = makeSineWaveImage(spatialFrequencyCPD, orientationDegs, X, Y)
    sigma = max(X(:))/3;
    gaussianEnvelope = exp(-0.5*X.^2/sigma^2) .* exp(-0.5*Y.^2/sigma^2);
    sinwaveImage = sin(2*pi*(X*spatialFrequencyCPD*cosd(orientationDegs) + Y*spatialFrequencyCPD*sind(orientationDegs)));
    sinwaveImage = sinwaveImage .* gaussianEnvelope;
end

function [psfImage, coneMosaicImage, psfSupportArcMin, Xhires,Yhires, NyquistCyclesPerDegreeFromConeSeparation] = ...
             makePSFandConeImage(cm, psf, coneApertureFraction, fftSize)
      
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
     
     % Interpolate to 1024 samples
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
    
    if (coneApertureFraction < 0.1)
        employDeltaFunctionsConeApertures = true;
    else
        employDeltaFunctionsConeApertures = ~true;
    end
    
    
    for iCone = 1:conesNum

        diffPos = bsxfun(@minus, XY,conePositionsMicrons(iCone,:));
        [~,idxCone] = min(sum(diffPos.^2,2));
            
        if (employDeltaFunctionsConeApertures)
            useActualConeApertures = ~true;
            idx = idxCone;
        else
            useActualConeApertures = true;
            xx = conePositionsMicrons(iCone,1) + coneApertureFraction*coneApertureMicrons(iCone)*0.5*coneOutlineX;
            yy = conePositionsMicrons(iCone,2) + coneApertureFraction*coneApertureMicrons(iCone)*0.5*coneOutlineY;
            [in,on] = inpolygon(Xhires(:), Yhires(:), xx,yy);
            idx = find((in|on));
        end
        
        % Add cone aperture
        tmp = Xhires*0;

        if (useActualConeApertures)
            tmp(idx) = 1;
        else
            tmp(idxCone) = 1;
        end
        coneMosaicImage = coneMosaicImage + tmp;
    end
    
end

function [coneMosaicActivationImage, stimulusImage, indicesList] = ...
             makeConeMosaicActivationImage(cm, coneApertureFraction, psfImage, stimulusSpatialFrequency, Xhires,Yhires, indicesList)
         
    t = clock;
         
    % Generate sinwave stimulus image
    orientationDegs = 0;
    stimulusImage = makeSineWaveImage(stimulusSpatialFrequency,orientationDegs, Xhires/cm.micronsPerDegree, Yhires/cm.micronsPerDegree);
    

    fftSize = size(stimulusImage,1);
    stimulusPSFImage = ifft2(fft2(stimulusImage, fftSize,fftSize) .* fft2(psfImage, fftSize,fftSize), fftSize,fftSize);
    
    % Generate sinewave stimulus optical image (convolution with PSF)
    %stimulusPSFImage = conv2(stimulusImage, psfImage, 'same');
    fprintf('Convolution in frequency domain took %2.1f seconds\n', etime(clock,t));

    % Centered cone positions
    conePositionsMicrons = bsxfun(@minus, cm.coneRFpositionsMicrons, cm.eccentricityMicrons);
    conesNum = size(conePositionsMicrons,1);
    % Cone apertures
    coneApertureMicrons = cMosaic.coneApertureToDiameterRatio * cm.coneRFspacingsMicrons;

    % Generate cone mosaic image
    coneOutlineX = cosd(0:15:360);
    coneOutlineY = sind(0:15:360);

    coneMosaicActivationImage = Xhires*0;
    
    XY = [Xhires(:) Yhires(:)];
    
    if (coneApertureFraction < 0.1)
        employDeltaFunctionsConeApertures = true;
    else
        employDeltaFunctionsConeApertures = ~true;
    end
    
    indicesListWasEmpty = isempty(indicesList);
    for iCone = 1:conesNum
        
        diffPos = bsxfun(@minus, XY,conePositionsMicrons(iCone,:));
        [~,idxCone] = min(sum(diffPos.^2,2));
        
        if (indicesListWasEmpty)
            if (employDeltaFunctionsConeApertures)
                %Use delta function cone apertures
                idx = idxCone;
            else
                % Use actual cone apertures
                xx = conePositionsMicrons(iCone,1) + coneApertureFraction*coneApertureMicrons(iCone)*0.5*coneOutlineX;
                yy = conePositionsMicrons(iCone,2) + coneApertureFraction*coneApertureMicrons(iCone)*0.5*coneOutlineY;
                [in,on] = inpolygon(Xhires(:), Yhires(:), xx,yy);
                idx = find((in|on));
            end
            indicesList{iCone} = idx;
        end
        
        
        % Add activation of cone computed as the spatial integration of the stimulus over the cone aperture
        tmp = Xhires*0;
        %tmp(idx) = sum(stimulusImage(idx));
        tmp(idxCone) = sum(stimulusPSFImage(indicesList{iCone}));
        coneMosaicActivationImage = coneMosaicActivationImage + tmp;
        
    end

    fprintf('Mosaic activation image took %2.1f seconds\n', etime(clock,t));
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

 
% Function to generate a semitransparent controur plot
function semiTransparentContourPlot(axesHandle, xSupport, ySupport, zData, zLevels, cmap, alpha, contourLineColor)
    % Compute contours at desired Z-level
    C = contourc(xSupport, ySupport, zData, zLevels);
    % Go through the contour matrix and plot each contour separately
    dataPoints = size(C,2);
    startPoint = 1;
    hold(axesHandle, 'on');
    while (startPoint < dataPoints)
        theLevel = C(1,startPoint);
        theNormalizedLevel = (theLevel-min(zLevels))/(max(zLevels)-min(zLevels));
        theLevelVerticesNum = C(2,startPoint);
        x = C(1,startPoint+(1:theLevelVerticesNum));
        y = C(2,startPoint+(1:theLevelVerticesNum));
        v = [x(:) y(:)];
        f = 1:numel(x);

        lutIndex = 1+round(theNormalizedLevel*(size(cmap,1)-1));
        patch('Faces', f, 'Vertices', v, 'EdgeColor', [0 0 0], ...
            'FaceColor', cmap(lutIndex,:), ...
            'FaceAlpha', alpha, ...
            'EdgeColor', contourLineColor, ...
            'LineStyle', '-', 'LineWidth', 1.0, ...
        'Parent', axesHandle);
        startPoint = startPoint + theLevelVerticesNum+1;
    end
end
 