function analyzePolansOptics(reAnalyzeData)

    % Get directory
    [directory,~] = fileparts(which(mfilename()));
    exportsDir = fullfile(directory, 'exports');
    
    if (reAnalyzeData)
        reAnalyze(exportsDir);
         % Rank subjects  
         rankStrategy = 'resolution'; %  Choose from {'resolution', 'peak resolution', 'correlation coefficient 10 degs'}
         rankedSubjectIDs = rankSubjects(exportsDir, rankStrategy)
         % Plot ranked subject data
         plotRankedSubjects(exportsDir, rankedSubjectIDs, rankStrategy);
        
    else
        doRankAnalysis = true;
        if (doRankAnalysis)
            % Plot ranked subject data
            rankStrategy = 'resolution'; %  Choose from {'resolution', 'peak resolution', 'correlation coefficient 10 degs'}
            rankedSubjectIDs = rankSubjects(exportsDir, rankStrategy)
            plotRankedSubjects(exportsDir, rankedSubjectIDs, rankStrategy);
        end
        
        plotSummary(exportsDir);
    end
end

function plotRankedSubjects(exportsDir, rankedSubjectIDs, rankStrategy)
    load(fullfile(exportsDir,'PolansOpticsAnalysis.mat'), 'subjectPSFData', 'coneCutoffSF', 'mosaicNyquistFrequencyCPD');
    
    subjectIDs = [];
    for includedSubjectCount = 1:numel(subjectPSFData)
        d = subjectPSFData{includedSubjectCount};
        horizontalEcc = d.horizontalEcc;
        verticalEcc = d.verticalEcc;
        subjectID = d.subjectID;
        psfXCutoffSF(includedSubjectCount,:,:) = d.psfXCutoffSF;
        psfYCutoffSF(includedSubjectCount,:,:) = d.psfYCutoffSF;
        subjectIDs = cat(2, subjectIDs, d.subjectID);
    end
    
    fovealEccY = find(verticalEcc == 0);
    psfXCutoffSF = squeeze(psfXCutoffSF(:,fovealEccY,:));
    psfYCutoffSF = squeeze(psfYCutoffSF(:,fovealEccY,:));
    meanPSFXCutoffSF = mean(psfXCutoffSF,1);
    
    rowsNum = 3;
    colsNum = 4;
    sv = NicePlot.getSubPlotPosVectors(...
           'colsNum', colsNum, ...
           'rowsNum', rowsNum, ...
           'heightMargin',  0.04, ...
           'widthMargin',    0.01, ...
           'leftMargin',     0.02, ...
           'rightMargin',    0.00, ...
           'bottomMargin',   0.04, ...
           'topMargin',      0.01); 
       
     hFig = figure();
     clf;
     set(hFig, 'Position', [10 10 2000 1400], 'Color', [1 1 1]);
     
     for includedSubjectCount = 1:numel(rankedSubjectIDs)
         
        subjectID = rankedSubjectIDs(includedSubjectCount);
        index = find(subjectIDs == rankedSubjectIDs(includedSubjectCount));
        
        r = floor((includedSubjectCount-1)/colsNum);
        r = mod(r,rowsNum)+1;
        c = mod(includedSubjectCount-1,colsNum)+1;

        subplot('Position', sv(r,c).v);
        
        plot(horizontalEcc, mosaicNyquistFrequencyCPD(fovealEccY,:), 'k-', 'LineWidth', 3.0); hold on;
        plot(horizontalEcc, mosaicNyquistFrequencyCPD(fovealEccY,:), 'g--', 'LineWidth', 1.0);
        if (strcmp(rankStrategy, 'correlation coefficient 10 degs'))
            plot(horizontalEcc, meanPSFXCutoffSF, 'm-', 'LineWidth', 1.5);
        end
        plot(horizontalEcc, psfXCutoffSF(index,:), 'bo-', 'MarkerFaceColor', [0.5 0.5 1], 'LineWidth', 1.5);
        plot(horizontalEcc, psfYCutoffSF(index,:), 'ro-', 'MarkerFaceColor', [1 0.5 0.5], 'LineWidth', 1.5);
        title(sprintf('subjectID: %d, rank: %d/%d', subjectID, includedSubjectCount,numel(rankedSubjectIDs)));
        
        axis 'square';
        grid on;
        set(gca, 'FontSize', 14, 'XLim', [-21 21], 'YLim', [0 70], 'XTick', -20:5:20, 'YTick', 0:5:70);
        
     end
     
     figTitle = sprintf('Ranked_POLANS');
     NicePlot.exportFigToPDF(fullfile(exportsDir,sprintf('%s.pdf', figTitle)), hFig, 300);
end


function rankedSubjectIDs = rankSubjects(exportsDir, rankStrategy)

    load(fullfile(exportsDir,'PolansOpticsAnalysis.mat'), 'subjectPSFData', 'coneCutoffSF', 'mosaicNyquistFrequencyCPD');
    
    subjectIDs = [];
    for includedSubjectCount = 1:numel(subjectPSFData)
        d = subjectPSFData{includedSubjectCount};
        horizontalEcc = d.horizontalEcc;
        verticalEcc = d.verticalEcc;
        subjectID = d.subjectID;
        psfXCutoffSF(includedSubjectCount,:,:) = d.psfXCutoffSF;
        psfYCutoffSF(includedSubjectCount,:,:) = d.psfYCutoffSF;
        subjectIDs = cat(2, subjectIDs, d.subjectID);
    end

    
    
    switch (rankStrategy)
        case 'correlation coefficient 10 degs'
            % Rank according to correlation coeff in [-10 10]
            idx = find(abs(horizontalEcc)<=10);
            fovealEccY = find(verticalEcc == 0);
            psfXCutoffSF = squeeze(psfXCutoffSF(:,fovealEccY,idx));
            psfYCutoffSF = squeeze(psfYCutoffSF(:,fovealEccY,idx));
            meanPSFXCutoffSF = mean(psfXCutoffSF,1);
            meanPSFYCutoffSF = mean(psfYCutoffSF,1);
            r = corr(meanPSFXCutoffSF', psfXCutoffSF');
            
        case 'peak resolution'
            fovealEccX = find(horizontalEcc == 0);
            fovealEccY = find(verticalEcc == 0);
            r = psfXCutoffSF(:,fovealEccY,fovealEccX);
            
        case 'resolution'
            fovealEccX = find(horizontalEcc == 0);
            fovealEccY = find(verticalEcc == 0);
            r = sqrt(psfXCutoffSF(:,fovealEccY,fovealEccX) .* psfYCutoffSF(:,fovealEccY,fovealEccX));
            
        otherwise
            error('Unknown rankStrategy')
    end
    
    [r,sortedSubjectIndices] = sort(r, 'descend');
    rankedSubjectIDs = subjectIDs(sortedSubjectIndices);
    
    hFig = figure();
    set(hFig, 'Color', [1 1 1], 'Position', [10 10 1200 800]);
    plot(1:numel(rankedSubjectIDs), r, 'bo-', 'MarkerFaceColor', [0.5 0.5 1], 'MarkerSize', 12, 'LineWidth', 1.5);
    xlabel(sprintf('subject ID'));
    switch (rankStrategy)
        case 'correlation coefficient 10 degs'
            ylabel('correlation coefficient');
            set(gca, 'YLim', [-1 1]);
        case 'peak resolution'
            ylabel('foveal resolution (c/deg)');
            set(gca, 'YLim', [0 65]);
    end
    set(gca, 'FontSize', 16, 'XLim',[0 numel(rankedSubjectIDs)+1], 'XTick', 1:numel(rankedSubjectIDs), 'XTickLabel', rankedSubjectIDs); 
    grid on
    NicePlot.exportFigToPDF(fullfile(exportsDir, 'PolansSubjectsRanked'), hFig, 300);
       
end


function plotSummary(exportsDir)
    load(fullfile(exportsDir,'PolansOpticsAnalysis.mat'), 'subjectPSFData', 'coneCutoffSF', 'mosaicNyquistFrequencyCPD');
    
    subjectsNum = numel(subjectPSFData);
    
    for includedSubjectCount = 1:subjectsNum
        d = subjectPSFData{includedSubjectCount};
        horizontalEcc = d.horizontalEcc;
        verticalEcc = d.verticalEcc;
        subjectID = d.subjectID;
        psfXCutoffSF(includedSubjectCount,:,:) = d.psfXCutoffSF;
        psfYCutoffSF(includedSubjectCount,:,:) = d.psfYCutoffSF;
        zCoeffs(includedSubjectCount,:,:,:) = d.zCoeffs;

        idxZeroVecc = find(verticalEcc == 0);
        [maxSF,idx] = max(squeeze(d.psfXCutoffSF(idxZeroVecc,:)));
        bestEcc(includedSubjectCount,:) = horizontalEcc(idx);
        sfCutoff(includedSubjectCount,:) = maxSF;
    end
    
    
    fovealEccX = find(horizontalEcc == 0);
    fovealEccY = find(verticalEcc == 0);
    
    
    rowsNum = 3;
    colsNum = 4;
    sv = NicePlot.getSubPlotPosVectors(...
           'colsNum', colsNum, ...
           'rowsNum', rowsNum, ...
           'heightMargin',  0.05, ...
           'widthMargin',    0.02, ...
           'leftMargin',     0.02, ...
           'rightMargin',    0.00, ...
           'bottomMargin',   0.04, ...
           'topMargin',      0.01); 
       
    theZs = squeeze(zCoeffs(:,fovealEccY,:,:));
    idx = find((horizontalEcc >= 13) & (horizontalEcc <= 18));
    theZs(:,idx,:) = nan;
  
    hFig = figure(111); clf;
    set(hFig, 'Position', [10 10 1400 975], 'Color', [1 1 1]);
    for zCoeffIndex = 4:15
        kk = zCoeffIndex-3;
        r = floor((kk-1)/colsNum);
        r = mod(r,rowsNum)+1;
        c = mod(kk-1,colsNum)+1;
        subplot('Position', sv(r,c).v);
        theZ2s = squeeze(theZs(:,:,zCoeffIndex));
        plot(horizontalEcc, theZ2s, 'bs-', 'LineWidth', 1.0);
        ylabel(sprintf('Z%d', zCoeffIndex-1));
        axis 'square';
        grid on;
        maxZamp = max([0.41 max(abs(theZ2s(:)))]);
        set(gca, 'FontSize', 12, 'YLim', maxZamp*[-1 1], 'XLim', [-21 21], 'XTick', -20:5:20);
        if (zCoeffIndex == 5)
            ylabel(sprintf('Z%d (defocus)', zCoeffIndex-1));
        elseif (zCoeffIndex == 4)
            ylabel(sprintf('Z%d (oblique astigmatism)', zCoeffIndex-1));
        elseif (zCoeffIndex == 6)
            ylabel(sprintf('Z%d (vertical astigmatism)', zCoeffIndex-1));
        end
        xlabel('eccentricity (degs)');
    end
    pause
    
    
    hFig = figure(1); clf;
    set(hFig, 'Position', [10 10 1300 510], 'Color', [1 1 1]);
    subplot(1,3,1);
    histogram(squeeze(zCoeffs(:,fovealEccY,fovealEccX, 4)), -1:0.1:1);
    xlabel('z3 (oblique astigmatism)');
    axis 'square'
    set(gca, 'XTick', -1:0.2:1, 'YLim', [0 6], 'YTick', [0:5]);
    
    subplot(1,3,2);
    histogram(squeeze(zCoeffs(:,fovealEccY,fovealEccX, 5)), -1:0.1:1);
    xlabel('z4 (defocus)');
    axis 'square'
    set(gca, 'XTick', -1:0.2:1, 'YLim', [0 6], 'YTick', [0:5]);
    
    subplot(1,3,3);
    histogram(squeeze(zCoeffs(:,fovealEccY,fovealEccX, 6)), -1:0.1:1);
    xlabel('z5 (vertical astigmatism)');
    axis 'square'
    set(gca, 'XTick', -1:0.2:1, 'YLim', [0 6], 'YTick', [0:5]);
    NicePlot.exportFigToPDF(fullfile(exportsDir,'PolansCoeffs.pdf'), hFig, 300);
    
    % Plot the position of best resolution
    hFig = figure(2); clf;
    set(hFig, 'Position', [10 10 950 410], 'Color', [1 1 1]);
    subplot(1,2,1);
    histogram(bestEcc, (-20:1:20)-0.5);
    text(14, 0.3, 'OD', 'FontSize', 12, 'FontWeight', 'Bold');
    xlabel('horizontal eccentricity of peak SF cutoff (degs)');
    ylabel('count');
    set(gca, 'FontSize', 12, 'XTick', -20:5:20, 'XLim', [-21 21]);
    axis 'square'
    
    subplot(1,2,2);
    plot(bestEcc, sfCutoff, 'k.');
    hold on;
    xlabel('horizontal eccentricity of peak SF cutoff(degs)');
    ylabel('peak SF cutoff (-3dB) (c/deg)');
    set(gca, 'FontSize', 12, 'XTick', -20:5:20, 'XLim', [-21 21], 'YLim', [0 65]);
    text(14, 4, 'OD', 'FontSize', 12, 'FontWeight', 'Bold');
    axis 'square'
    NicePlot.exportFigToPDF(fullfile(exportsDir,'PolansBestEcc.pdf'), hFig, 300);
    
    
    
    rowsNum = 3;
    colsNum = 4;
    sv = NicePlot.getSubPlotPosVectors(...
           'colsNum', colsNum, ...
           'rowsNum', rowsNum, ...
           'heightMargin',  0.05, ...
           'widthMargin',    0.01, ...
           'leftMargin',     0.02, ...
           'rightMargin',    0.00, ...
           'bottomMargin',   0.03, ...
           'topMargin',      0.02); 
       
    hFig = figure(100); clf
    set(hFig, 'Color', [1 1 1], 'Position', [10 10 2000 1400]);

    for vIndex = 1:size(mosaicNyquistFrequencyCPD,1)
            
            if (vIndex <=4)
                r = 1;
                c = vIndex;
            elseif (vIndex <= 7)
                r = 2;
                c = vIndex-4;
            else
                r = 3;
                c = vIndex-7;
            end
        subplot('Position', sv(r,c).v);
        
        h1 = plot(horizontalEcc, mosaicNyquistFrequencyCPD(vIndex,:), 'k-', 'LineWidth', 3);
        hold on
        plot(horizontalEcc, mosaicNyquistFrequencyCPD(vIndex,:), 'g--', 'LineWidth', 1.5);
        h3 = plot(horizontalEcc, mean(squeeze(psfXCutoffSF(:, vIndex,:)),1), 'bo-', 'LineWidth', 1.5, 'MarkerSize', 10, 'MarkerFaceColor', [0.5 0.5 1]);
        h4 = plot(horizontalEcc, mean(squeeze(psfYCutoffSF(:, vIndex,:)),1), 'ro-', 'LineWidth', 1.5, 'MarkerSize', 10, 'MarkerFaceColor', [1 0.5 0.5]);
        if (verticalEcc(vIndex) == 0)
        legend([h1 h3 h4], {'mosaic Nyquist freq.', 'PSF_x (mean)', 'PSF_y (mean)'});
        xlabel('horizontal eccentricity (deg)');
        end
        
        %ylabel('spatial frequency cutoff, -15dB (c/deg)');
        title(sprintf('%2.0f deg', verticalEcc(vIndex)));
        axis 'square';  grid 'on'; box 'off'
        set(gca, 'XLim', [-20 20], 'YLim', [0 70], 'XTick', -20:5:20, 'YTick', 0:5:100,  'FontSize', 16);
        set(gca, 'Color', 'none', 'XColor', 'none', 'YColor', 'none', 'LineWidth', 1.0);
        if (verticalEcc(vIndex) == 0)
            set(gca, 'Color', 'none', 'XColor', [0.3 0.3 0.3], 'YColor', [0.3 0.3 0.3]);
        end
        
        
    end
    NicePlot.exportFigToPDF(fullfile(exportsDir,sprintf('PolansAcrossAllEccentricities.pdf')), hFig, 300);
    
    
    vIndex = find(verticalEcc == 0);
    
    hFig = figure(2);
    clf;
    set(hFig, 'Color', [1 1 1], 'Position', [10 10 1500 710]);
    ax = subplot(1,2,1);
    h1 = plot(ax, horizontalEcc, mosaicNyquistFrequencyCPD(vIndex,:), 'k-', 'LineWidth', 3);
    hold(ax, 'on');
    h2 = plot(ax,horizontalEcc, mosaicNyquistFrequencyCPD(vIndex,:), 'g--', 'LineWidth', 1.5);
    h3 = plot(ax,horizontalEcc, squeeze(mean(psfXCutoffSF(:, vIndex,:),1)), 'b-', 'LineWidth', 1.5);
    h4 = shadeAreaBetweenCyrves(ax, horizontalEcc, squeeze(max(psfXCutoffSF(:, vIndex,:),[],1)), squeeze(min(psfXCutoffSF(:, vIndex,:),[],1)), [0.5 0.5 1], 0.5);
    legend([h1 h3], {'mosaic Nyquist freq.', 'PSF cutoff (mean)'});
    xlabel('horizontal eccentricity (deg)');
    ylabel('spatial frequency cutoff, -15dB (c/deg)');
    title('high-resolution axis');
    axis 'square';  grid 'on'; box 'off'
    set(gca, 'XLim', [-20 20], 'YLim', [0 70], 'XTick', -20:5:20, 'YTick', 0:5:100,  'FontSize', 16);
    set(gca, 'XColor', [0.3 0.3 0.3], 'YColor', [0.3 0.3 0.3], 'LineWidth', 1.0)
    
    ax = subplot(1,2,2);
    h1 = plot(ax,horizontalEcc, mosaicNyquistFrequencyCPD(vIndex,:), 'k-', 'LineWidth', 3);
    hold(ax, 'on');
    h2 = plot(ax,horizontalEcc, mosaicNyquistFrequencyCPD(vIndex,:), 'g--', 'LineWidth', 1.5);
    h3 = plot(ax,horizontalEcc, squeeze(mean(psfYCutoffSF(:, vIndex,:),1)), 'r-', 'LineWidth', 1.5);
    h4 = shadeAreaBetweenCyrves(ax, horizontalEcc, squeeze(max(psfYCutoffSF(:, vIndex,:),[],1)), squeeze(min(psfYCutoffSF(:, vIndex,:),[],1)), [1 0.5 0.5], 0.5);
    legend([h1 h3], {'mosaic Nyquist freq.', 'PSF cutoff (mean)'});
    xlabel('horizontal eccentricity (deg)');
    ylabel('spatial frequency cutoff, -15dB (c/deg)');
    title('low-resolution axis');
    axis 'square';  grid 'on'; box 'off'
    set(gca, 'XLim', [-20 20], 'YLim', [0 70], 'XTick', -20:5:20, 'YTick', 0:5:100, 'FontSize', 16);
    set(gca, 'XColor', [0.3 0.3 0.3], 'YColor', [0.3 0.3 0.3], 'LineWidth', 1.0)
    NicePlot.exportFigToPDF(fullfile(exportsDir,sprintf('PolansAcrossHorizontalEccentricity.pdf')), hFig, 300);
end

function h = shadeAreaBetweenCyrves(ax, x, curve1, curve2, shadeColor, alphaValue)
    x = x(:);
    curve1 = curve1(:);
    curve2 = curve2(:);
    x2 = [x', fliplr(x')];
    inBetween = [curve1', fliplr(curve2')];
    size(inBetween)
    size(x2)
    h = fill(ax,x2, inBetween, shadeColor, 'FaceAlpha', alphaValue);
end


function reAnalyze(exportsDir)
    
    opticsParams = struct(...
        'zernikeDataBase',  'Polans2015', ...
        'subjectID',1, ...
        'subtractCentralRefraction', false, ...
        'zeroCenterPSF', true, ...
        'flipPSFUpsideDown', true, ...
        'whichEye', 'right eye', ...
        'pupilDiameterMM', 3);
   
    
    horizontalEcc = -25:25;
    verticalEcc = PolansOptics.constants.measurementVerticalEccentricities;
    
    % Subject indices
    subjectIndices = 1:10;
    
    whichEyes = {'right eye'};

    plotEachPosition = ~true;
    
    rowsNum = 3;
    colsNum = 4;
    sv = NicePlot.getSubPlotPosVectors(...
           'colsNum', colsNum, ...
           'rowsNum', rowsNum, ...
           'heightMargin',  0.05, ...
           'widthMargin',    0.01, ...
           'leftMargin',     0.02, ...
           'rightMargin',    0.00, ...
           'bottomMargin',   0.03, ...
           'topMargin',      0.02); 
     
    if (plotEachPosition)
        videoFileName = fullfile(exportsDir,sprintf('%s_Subject%d_%s_Pupil%3.2fmm_SubtractCentralRefraction', ...
                    opticsParams.zernikeDataBase,  opticsParams.subjectID, opticsParams.whichEye, opticsParams.pupilDiameterMM));
        videoOBJ = VideoWriter(videoFileName, 'MPEG-4');
        videoOBJ.FrameRate = 30;
        videoOBJ.Quality = 100;
        videoOBJ.open();
    end
    
    
    eyeIndex = 1;
    
    for subjectID = 1:numel(subjectIndices)
        
        hFig = figure(subjectID);
        clf;
        set(hFig, 'Position', [10 10 2000 1400]);
     
        for vEcc = 1:numel(verticalEcc)
            
            fprintf('ecc:%2.1f degs, subject %d, subtract central refraction: %d\n', verticalEcc(vEcc), subjectIndices(subjectID), opticsParams.subtractCentralRefraction);
            
            if (vEcc <=4)
                r = 1;
                c = vEcc;
            elseif (vEcc <= 7)
                r = 2;
                c = vEcc-4;
            else
                r = 3;
                c = vEcc-7;
            end
            
            for iEcc = 1:numel(horizontalEcc) 

                % Update params
                mosaicEcc = [horizontalEcc(iEcc) verticalEcc(vEcc)];
                opticsParams.subjectID = subjectIndices(subjectID);
                opticsParams.whichEye = whichEyes{eyeIndex};
                opticsParams.subtractCentralRefraction = PolansOptics.constants.subjectRequiresCentralRefractionCorrection(opticsParams.subjectID);
                
                % Compute PSFimage and cone image (single aperture)
                [psfImage, coneMosaicImage, NyquistFrequency, psfSupportArcMin, theZCoeffs] = psfAndConeImage(mosaicEcc, opticsParams);
                
                if (isempty(psfImage))
                    fprintf('No data for ecc (x:%2.0f, y:%2.0f)\n', horizontalEcc(iEcc), verticalEcc(vEcc));
                    continue;
                end

                mosaicNyquistFrequencyCPD(vEcc, iEcc) = NyquistFrequency;
                zCoeffs(vEcc, iEcc,:) = theZCoeffs;
                
                % Compute power spectra
                [psfImageSpectrum, psfImageSpectrumRoatated, coneMosaicImageSpectrum, ...
                    spectralSupportCyclesPerDegree, spectralSupportCyclesPerDegreePositive, coneSpectrumSlice, ...
                    spectralSupportCyclesPerDegreePositiveX, psfSpectrumSliceX, ...
                    spectralSupportCyclesPerDegreePositiveY, psfSpectrumSliceY, ...
                    effectivePSFSpectrumSliceX, effectivePSFSpectrumSliceY, ...
                    coneCutoffSF(vEcc, iEcc), psfXCutoffSF(vEcc, iEcc), psfYCutoffSF(vEcc, iEcc)] = ...
                    computePSFAndConePowerSpectra(psfImage, coneMosaicImage, psfSupportArcMin);
        
                if (plotEachPosition)
                    hFigVideo = visualizeSpectralAnalysis(1000,horizontalEcc(iEcc), verticalEcc(vEcc), ...
                        psfSupportArcMin, psfImage, coneMosaicImage, ...
                        spectralSupportCyclesPerDegree, spectralSupportCyclesPerDegreePositive, ...
                        psfImageSpectrum, psfImageSpectrumRoatated,  ...
                        coneMosaicImageSpectrum, coneSpectrumSlice, ...
                        spectralSupportCyclesPerDegreePositiveX, psfSpectrumSliceX, effectivePSFSpectrumSliceX,...
                        spectralSupportCyclesPerDegreePositiveY, psfSpectrumSliceY, effectivePSFSpectrumSliceY,...
                        coneCutoffSF(vEcc, iEcc), psfXCutoffSF(vEcc, iEcc), psfYCutoffSF(vEcc, iEcc), ...
                        mosaicNyquistFrequencyCPD(vEcc, iEcc));
                    % Add frame to video
                    drawnow;
                    videoOBJ.writeVideo(getframe(hFigVideo)); 
                end
           
                hFig = figure(subjectID);
                subplot('Position', sv(r,c).v);


                plot(horizontalEcc(1:iEcc), mosaicNyquistFrequencyCPD(vEcc, 1:iEcc), 'ko-', 'LineWidth', 1.5, 'MarkerSize', 8, 'MarkerFaceColor', [0.8 0.8 0.8]); 
                hold on;
                plot(horizontalEcc(1:iEcc), psfXCutoffSF(vEcc, 1:iEcc), 'bo-', 'LineWidth', 1.5, 'MarkerSize', 8, 'MarkerFaceColor', [0.5 0.5 1]);
                plot(horizontalEcc(1:iEcc), psfYCutoffSF(vEcc, 1:iEcc), 'ro-', 'LineWidth', 1.5, 'MarkerSize', 8, 'MarkerFaceColor', [1 0.5 0.5]);
                xlabel('horizontal eccentricity (deg)');
                ylabel('spatial frequency (c/deg)');
                title(sprintf('vertical eccentricity: %2.0f', verticalEcc(vEcc)));
                legend({'mosaic Nyquist freq.', 'psf_x cutoff (-15dB)', 'psf_y cutoff (-15dB)'});

                axis 'square';
                grid on;
                set(gca, 'FontSize', 14, 'XLim', [-26 26], 'YLim', [0 70], 'XTick', -20:5:20, 'YTick', 0:5:70);

            end
             
            fprintf('Finished with subject %d\n', opticsParams.subjectID);
        end
        
        
        % Data for exporting
        subjectPSFData{subjectID} = struct(...
            'horizontalEcc', horizontalEcc, ...
            'verticalEcc', verticalEcc, ...
            'subjectID', opticsParams.subjectID, ...
            'zCoeffs', zCoeffs, ...
            'whichEye', opticsParams.whichEye, ...
            'pupilDiamMM', opticsParams.pupilDiameterMM, ...
            'psfXCutoffSF', psfXCutoffSF, ...
            'psfYCutoffSF', psfYCutoffSF ...
            );
                
        if (plotEachPosition)
            videoOBJ.close();
        end
        
        if (opticsParams.subtractCentralRefraction)
            figTitle = sprintf('%s_Subject%d_%s_Pupil%3.2fmm_SubtractCentralRefraction', ...
                opticsParams.zernikeDataBase,  opticsParams.subjectID, opticsParams.whichEye, opticsParams.pupilDiameterMM);
        else
            figTitle = sprintf('%s_Subject%d_%s_Pupil%3.2fmm_DoNotSubtractCentralRefraction', ...
                opticsParams.zernikeDataBase, opticsParams.subjectID, opticsParams.whichEye, opticsParams.pupilDiameterMM);
        end
        NicePlot.exportFigToPDF(fullfile(exportsDir,sprintf('%s.pdf', figTitle)), hFig, 300);
    end
    
    
    % Export data
    save(fullfile(exportsDir,'PolansOpticsAnalysis.mat'), 'subjectPSFData', 'coneCutoffSF', 'mosaicNyquistFrequencyCPD');
    
end


