function analyzeOptics()

   
    opticsParams = struct(...
        'zernikeDataBase',  'Polans2015', ...
        'subjectID',10, ...
        'subtractCentralRefraction', true, ...
        'zeroCenterPSF', ~true, ...
        'flipPSFUpsideDown', true, ...
        'whichEye', 'right eye', ...
        'pupilDiameterMM', 3);
    
   %%% s#86613002
    
    
    horizontalEcc = -8:0; % :25;
    verticalEcc = horizontalEcc*0 + 0;
    whichEyes = {'right eye'}; % {'left eye', 'right eye'};
    
    plotEachPosition = true;
    
    for eyeIndex = 1:numel(whichEyes)
    for subjectID = 7:7 % 1:10
    for subtractCentralRefraction = 0:0
        for iEcc = 1:numel(horizontalEcc)
 
            fprintf('ecc:%2.1f degs, subject %d, subtract central refraction: %d\n', horizontalEcc(iEcc), subjectID, subtractCentralRefraction);
            
            % Update params
            mosaicEcc = [horizontalEcc(iEcc) verticalEcc(iEcc)];
            opticsParams.subjectID = subjectID;
            opticsParams.whichEye = whichEyes{eyeIndex};
            opticsParams.subtractCentralRefraction = (subtractCentralRefraction==1);
            
            % Compute PSFimage and cone mosaic image
            [psfImage, coneMosaicImage, psfSupportArcMin] = psfAndConeImage(mosaicEcc, opticsParams);
            

            % Compute power spectra
            [psfImageSpectrum, psfImageSpectrumRoatated, coneMosaicImageSpectrum, ...
                spectralSupportCyclesPerDegree, spectralSupportCyclesPerDegreePositive, coneSpectrumSlice, ...
                spectralSupportCyclesPerDegreePositiveX, psfSpectrumSliceX, ...
                spectralSupportCyclesPerDegreePositiveY, psfSpectrumSliceY, ...
                coneCutoffSF(iEcc), psfXCutoffSF(iEcc), psfYCutoffSF(iEcc)] = ...
                computePowerSpectra(psfImage, coneMosaicImage, psfSupportArcMin);

            coneCutoffSF(iEcc)
            
            if (plotEachPosition)
                visualizeSpectralAnalysis(iEcc,horizontalEcc(iEcc), verticalEcc(iEcc), ...
                    psfSupportArcMin, psfImage, coneMosaicImage, ...
                    spectralSupportCyclesPerDegree, spectralSupportCyclesPerDegreePositive, ...
                    psfImageSpectrum, psfImageSpectrumRoatated,  ...
                    coneMosaicImageSpectrum, coneSpectrumSlice, ...
                    spectralSupportCyclesPerDegreePositiveX, psfSpectrumSliceX, ...
                    spectralSupportCyclesPerDegreePositiveY, psfSpectrumSliceY, ...
                    coneCutoffSF(iEcc), psfXCutoffSF(iEcc), psfYCutoffSF(iEcc));
            end
        end

        hFig = figure();
        clf;
        set(hFig, 'Position', [10 10 1100 800]);
        plot(horizontalEcc, coneCutoffSF, 'ko-', 'LineWidth', 1.5, 'MarkerSize', 12, 'MarkerFaceColor', [0.8 0.8 0.8]); hold on;
        plot(horizontalEcc, psfXCutoffSF, 'bo-', 'LineWidth', 1.5, 'MarkerSize', 12, 'MarkerFaceColor', [0.5 0.5 1]);
        plot(horizontalEcc, psfYCutoffSF, 'ro-', 'LineWidth', 1.5, 'MarkerSize', 12, 'MarkerFaceColor', [1 0.5 0.5]);
        xlabel('horizontal eccentricity (deg)');
        ylabel('spatial frequency cutoff (-3dB), c/deg');
        legend({'cone','psf_x', 'psf_y'});
        
        axis 'square';
        grid on;
        set(gca, 'FontSize', 14, 'XLim', [-26 26], 'YLim', [0 30], 'XTick', -20:5:20, 'YTick', 0:5:50);

        if (opticsParams.subtractCentralRefraction)
            figTitle = sprintf('%s_Subject%d_%s_Pupil%3.2fmm_SubtractCentralRefraction', ...
                opticsParams.zernikeDataBase,  opticsParams.subjectID, opticsParams.whichEye, opticsParams.pupilDiameterMM);
        else
            figTitle = sprintf('%s_Subject%d_%s_Pupil%3.2fmm_DoNotSubtractCentralRefraction', ...
                opticsParams.zernikeDataBase, opticsParams.subjectID, opticsParams.whichEye, opticsParams.pupilDiameterMM);
        end
        title(figTitle);
        NicePlot.exportFigToPDF(sprintf('%s.pdf', figTitle), hFig, 300);
    end
    end
    end
    
end


function [psfImageSpectrum, psfImageSpectrumRoatated, coneMosaicImageSpectrum, ...
    spectralSupportCyclesPerDegree, spectralSupportCyclesPerDegreePositive, coneSpectrumSlice, ...
    spectralSupportCyclesPerDegreePositiveX, psfSpectrumSliceX, ...
    spectralSupportCyclesPerDegreePositiveY, psfSpectrumSliceY, ...
    coneCutoffSF, psfXCutoffSF, psfYCutoffSF] = ...
        computePowerSpectra(psfImage, coneMosaicImage, psfSupport)
    
    fftSize = 4096;
    psfImageSpectrum = fftshift(abs(fft2(psfImage, fftSize,fftSize)));
    psfImageSpectrum = psfImageSpectrum / max(psfImageSpectrum(:));
    
    binaryImage = psfImageSpectrum;
    m1 = min(binaryImage(:));
    m2 = max(binaryImage(:));
    binaryImage = imbinarize((binaryImage - m1)/(m2-m1));
    s = regionprops(binaryImage,'orientation');
    psfRotation = s.Orientation;
    
    psfImageSpectrumRoatated = imrotate(psfImageSpectrum, -psfRotation, 'bilinear', 'crop');
    psfImageSpectrumRoatated = psfImageSpectrumRoatated / max(psfImageSpectrumRoatated(:));
    
    
    
    dx = psfSupport(2)-psfSupport(1);
    maxSF = 1/(2*dx);
    deltaSF = maxSF / (fftSize/2);
    spectralSupportCyclesPerDegree = ((-maxSF+deltaSF):deltaSF:maxSF) * 60;
    
    N = size(psfImageSpectrumRoatated,1);
    midPoint = N/2+1;
    spectralSupportCyclesPerDegreePositive = spectralSupportCyclesPerDegree(midPoint(1):end);
    
    
    [maxPSF,maxPoint] = max(psfImageSpectrumRoatated(:),[], 1);
    [maxPoint(1), maxPoint(2)] = ind2sub(size(psfImageSpectrumRoatated), maxPoint);
    
    if (midPoint ~= maxPoint(1))
        fprintf('midPoint (%d) differs from maxPointY (%d)\n', midPoint,maxPoint(1));
    end
    if (midPoint ~= maxPoint(2))
        fprintf('midPoint (%d) differs from maxPointX (%d)\n', midPoint,maxPoint(2));
    end
    midPoint = maxPoint;
    
    
    spectralSupportCyclesPerDegreePositiveY = spectralSupportCyclesPerDegree(midPoint(1):end);
    psfSpectrumSliceY = psfImageSpectrumRoatated(midPoint(1):end,midPoint(2));
    
    spectralSupportCyclesPerDegreePositiveX = spectralSupportCyclesPerDegree(midPoint(2):end);
    psfSpectrumSliceX = psfImageSpectrumRoatated(midPoint(1), midPoint(2):end);
    

    % -3.01 dB in amplitude is half power
    cornetFrequencyAttenuationDB = -3.01;
    cornerFrequencyAttenuation = 10^(cornetFrequencyAttenuationDB/20);
    
    
    [d, idx] = sort(abs(psfSpectrumSliceX-cornerFrequencyAttenuation));
    w(1) = d(2)/(d(1)+d(2));
    w(2) = d(1)/(d(1)+d(2));
    psfXCutoffSF = spectralSupportCyclesPerDegreePositiveX(idx(1))*w(1) + spectralSupportCyclesPerDegreePositiveX(idx(2))*w(2);
    [d, idx] = sort(abs(psfSpectrumSliceY-cornerFrequencyAttenuation));
    w(1) = d(2)/(d(1)+d(2));
    w(2) = d(1)/(d(1)+d(2));
    psfYCutoffSF = spectralSupportCyclesPerDegreePositiveY(idx(1))*w(1) + spectralSupportCyclesPerDegreePositiveY(idx(2))*w(2);
    
    if (~isempty(coneMosaicImage))
        deltaAngle = 360;
        for k = 0:deltaAngle:(360-deltaAngle)
            rotatedConeMosaicImage = imrotate(coneMosaicImage, k, 'bilinear', 'crop');
            if (k == 0)
                coneMosaicImageSpectrum = fftshift(abs(fft2(rotatedConeMosaicImage, fftSize,fftSize)));
            else
                coneMosaicImageSpectrum = coneMosaicImageSpectrum + fftshift(abs(fft2(rotatedConeMosaicImage, fftSize,fftSize)));
            end
        end
        coneMosaicImageSpectrum = coneMosaicImageSpectrum / max(coneMosaicImageSpectrum(:));
        midPoint = N/2+1;
        coneSpectrumSlice = coneMosaicImageSpectrum(midPoint,midPoint:end);
        [d,idx] = sort(abs(coneSpectrumSlice-cornerFrequencyAttenuation));
        w(1) = d(2)/(d(1)+d(2));
        w(2) = d(1)/(d(1)+d(2));
        coneCutoffSF = spectralSupportCyclesPerDegreePositive(idx(1))*w(1) + spectralSupportCyclesPerDegreePositive(idx(2))*w(2);
    else
        coneMosaicImageSpectrum = [];  
        coneSpectrumSlice = nan(size(spectralSupportCyclesPerDegreePositive));  
        coneCutoffSF = nan;
    end
    
end


function [psfImage, coneMosaicImage, psfSupportArcMin] = psfAndConeImage(mosaicEcc, opticsParams)

    mosaicSizeDegs = max([0.1 0.15*sqrt(sum(mosaicEcc.^2,2))]);
    % Generate mosaic centered at target eccentricity
    cm = cMosaic(...
        'sizeDegs', mosaicSizeDegs*[1 1], ...    % SIZE in degs
        'eccentricityDegs', mosaicEcc, ...  % ECC in degs
        'opticalImagePositionDegs', 'mosaic-centered', ...
        'whichEye', opticsParams.whichEye ...
        );
    
    wavelength = 550;
    [oiEnsemble, psfEnsemble] = ...
            cm.oiEnsembleGenerate(mosaicEcc, ...
            'zernikeDataBase', opticsParams.zernikeDataBase, ...
            'subjectID', opticsParams.subjectID, ...
            'pupilDiameterMM', opticsParams.pupilDiameterMM, ...
            'subtractCentralRefraction', opticsParams.subtractCentralRefraction, ...
            'wavefrontSpatialSamples', 701, ...
            'zeroCenterPSF', opticsParams.zeroCenterPSF, ...
            'flipPSFUpsideDown', opticsParams.flipPSFUpsideDown);
        
    % Make PSF image
    psf = psfEnsemble{1};
    [~,idx] = min(abs(psf.supportWavelength-wavelength));
    psfImage = squeeze(psf.data(:,:,idx));
     
    psfSupportMicrons = cm.micronsPerDegree * psf.supportX/60;
    [Xorig,Yorig] = meshgrid(psfSupportMicrons,psfSupportMicrons);
     
    interpolationSamples = 3000;
    psfSupportMicrons = linspace(psfSupportMicrons(1), psfSupportMicrons(end), interpolationSamples);
    [Xhires,Yhires] = meshgrid(psfSupportMicrons,psfSupportMicrons);
    psfImage = interp2(Xorig, Yorig, psfImage, Xhires, Yhires);
    psfImage = psfImage / max(psfImage(:));
    psfSupportArcMin = psfSupportMicrons / cm.micronsPerDegree * 60;
    
    % Make cone mosaic image
    conePositions = cm.coneRFpositionsMicrons;
    if (isempty(cm.coneRFspacingsMicrons))
        coneMosaicImage = [];
        return;
    end
    
    conePositions = bsxfun(@minus, conePositions, mean(conePositions,1));
    coneRadialEcc = sqrt(sum(conePositions.^2,2));
    [~,idx] = sort(coneRadialEcc);

    coneRFspacing = mean(cm.coneRFspacingsMicrons(idx(1:25)));
    coneRFDiameter = cMosaic.coneApertureToDiameterRatio * coneRFspacing;

    % Just the center cone and its 6 neighbors
    for k = 0:6
        conePosition = [cosd(k*60+30) sind(k*60+30)] * coneRFspacing;
        if (k == 0)
            xx = 0;
            yy = 0;
        else
            xx = conePosition(1);
            yy = conePosition(2);
        end
        xx = xx + coneRFDiameter/2 * cosd(0:5:360);
        yy = yy + coneRFDiameter/2 * sind(0:5:360);
        [in,on] = inpolygon(Xhires(:), Yhires(:), xx,yy);
        idx = find((in|on));
        tmp = Xhires*0;
        tmp(idx) = 1;
        if (k == 0)
            coneMosaicImage = tmp;
        else
            coneMosaicImage = coneMosaicImage + tmp;
        end
     end
     coneMosaicImage = coneMosaicImage / max(coneMosaicImage(:));
     
     
end

function visualizeSpectralAnalysis(figNo,  horizontalEcc, verticalEcc, ...
    psfSupportArcMin, psfImage, coneMosaicImage, ...
    spectralSupportCyclesPerDegree, spectralSupportCyclesPerDegreePositive, ...
    psfImageSpectrum, psfImageSpectrumRoatated,  ...
    coneMosaicImageSpectrum, coneSpectrumSlice, ...
    spectralSupportCyclesPerDegreePositiveX, psfSpectrumSliceX, ...
    spectralSupportCyclesPerDegreePositiveY, psfSpectrumSliceY, ...
    coneCutoffSF, psfXCutoffSF, psfYCutoffSF)

    subplotPosVectors = NicePlot.getSubPlotPosVectors(...
        'rowsNum', 3, ...
        'colsNum', 2, ...
        'heightMargin',  0.06, ...
        'widthMargin',    0.03, ...
        'leftMargin',     0.04, ...
        'rightMargin',    0.01, ...
        'bottomMargin',   0.04, ...
        'topMargin',      0.01);

    hFig = figure(figNo); clf;
    set(hFig, 'Position', [10 10 870 1300]);

    % The PSF
    ax = subplot('Position', subplotPosVectors(1,1).v);
    imagesc(ax,psfSupportArcMin, psfSupportArcMin, psfImage);
    axis(ax,'image')
    set(ax, 'XLim', 20*[-1 1], 'YLim', 20*[-1 1], 'XTick', -20:10:20, 'YTick', -20:10:20, 'FontSize', 14);
    ylabel('space,y (arc min)');
    xlabel('space, x (arc min)');
    title(ax, sprintf('PSF @ (%2.0f,%2.0f)', horizontalEcc, verticalEcc));
    colormap(ax, brewermap(1024, '*spectral'));

    if (~isempty(coneMosaicImage))
        %The cone mosaic
        ax = subplot('Position', subplotPosVectors(1,2).v);
        imagesc(ax, psfSupportArcMin, psfSupportArcMin, coneMosaicImage);
        axis(ax, 'image')
        xlabel('space, x(arc min)');
        set(ax, 'XLim', 20*[-1 1], 'YLim', 20*[-1 1], 'XTick', -20:10:20, 'YTick', -20:10:20, 'FontSize', 14);
        title(ax, 'cone mosaic');
        colormap(ax, brewermap(1024, '*greys'));
    end
    
    % The amplitude spectrum of the PSF
    ax = subplot('Position', subplotPosVectors(2,1).v);
    imagesc(ax,spectralSupportCyclesPerDegree, spectralSupportCyclesPerDegree, log10(psfImageSpectrum.^2));
    axis(ax, 'image');
    set(ax, 'XLim', [-100 100], 'YLim', [-100 100], 'XTick', -80:20:80, 'YTick', -80:20:80, 'CLim', [-3 0], 'FontSize', 14);
    xlabel(ax,'spatial frequency, x (cpd)');
    ylabel(ax,'spatial frequency, y (cpd)');
    title(ax,'PSF power spectrum');
    colormap(ax,brewermap(1024, '*greys'));

    if (~isempty(coneMosaicImageSpectrum))
        % The amplitude spectrum of the cone mosaic
        ax = subplot('Position', subplotPosVectors(2,2).v);
        imagesc(ax,spectralSupportCyclesPerDegree, spectralSupportCyclesPerDegree, log10(coneMosaicImageSpectrum.^2));
        axis(ax, 'image')
        xlabel(ax,'spatial frequency, x (cpd)');
        set(ax, 'XLim', [-100 100], 'YLim', [-100 100], 'XTick', -80:20:80, 'YTick', -80:20:80, 'CLim', [-3 0],'FontSize', 14);
        title(ax,'cone mosaic power spectrum');
        colormap(ax,brewermap(1024, '*greys'));
    end
    
    % The amplitude spectrum of the rotated PSF
    ax = subplot('Position', subplotPosVectors(3,1).v);
    imagesc(ax,spectralSupportCyclesPerDegree, spectralSupportCyclesPerDegree, log10(psfImageSpectrumRoatated.^2));
    axis(ax, 'image')
    set(ax, 'XLim', [-100 100], 'YLim', [-100 100], 'XTick', -80:20:80, 'YTick', -80:20:80, 'CLim', [-3 0],'FontSize', 14);
    xlabel(ax,'spatial frequency, x (cpd)');
    ylabel(ax,'spatial frequency, y (cpd)');
    title(ax,'rotated PSF power spectrum');
    colormap(ax,brewermap(1024, '*greys'));

    
    % The slices of the power spectra
    ax = subplot('Position', subplotPosVectors(3,2).v);
    plot(ax, spectralSupportCyclesPerDegreePositive, coneSpectrumSlice.^2, 'ko-', ...
        'MarkerSize', 8, 'MarkerFaceColor', [0.8 0.8 0.8], 'LineWidth', 1.5); hold on;
    if (spectralSupportCyclesPerDegreePositiveX(1) < 0.3)
        spectralSupportCyclesPerDegreePositiveX(1) = 0.3;
    end
    if (spectralSupportCyclesPerDegreePositiveY(1) < 0.3)
        spectralSupportCyclesPerDegreePositiveY(1) = 0.3;
    end
    plot(ax, spectralSupportCyclesPerDegreePositiveX, psfSpectrumSliceX.^2, 'bo-', ...
        'MarkerSize', 8, 'MarkerFaceColor', [0.5 0.5 1], 'LineWidth', 1.5);
    plot(ax, spectralSupportCyclesPerDegreePositiveY, psfSpectrumSliceY.^2, 'ro-', ...
        'MarkerSize', 8, 'MarkerFaceColor', [1 0.5 0.5], 'LineWidth', 1.5);

    scatter(ax, coneCutoffSF, 0.5, ...
        14*14, 'MarkerEdgeColor', [ 0 0 0], 'MarkerFaceColor', [0.8 0.8 0.8], 'MarkerFaceAlpha', 0.5);
    scatter(ax, psfXCutoffSF, 0.5, ...
        14*14, 'MarkerEdgeColor', [ 0 0 1], 'MarkerFaceColor', [0.5 0.5 1], 'MarkerFaceAlpha', 0.5);
    scatter(ax, psfYCutoffSF, 0.5, ...
        14*14, 'MarkerEdgeColor', [1 0 0], 'MarkerFaceColor', [1 0.5 0.5], 'MarkerFaceAlpha', 0.5);
    set(ax, 'XLim', [0.3 100], 'XTick', [0.1 0.3 1 3 10 30 100], ...
        'YLim', [0 1], 'XScale', 'log', 'YScale', 'linear', 'FontSize', 14);
    legend(ax,{'cone', 'psf_x', 'psf_y'});
    title(ax,'power spectra slices');
    axis(ax, 'square');
    grid(ax,'on');
    xlabel(ax,'spatial frequency, x (cpd)');
    ylabel(ax,'power spectrum (norm.)');


    drawnow;
end
