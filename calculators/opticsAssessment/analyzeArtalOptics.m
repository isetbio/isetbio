function analyzeArtalOptics(reAnalyzeData)
   
    % Get directory
    [directory,~] = fileparts(which(mfilename()));
    exportsDir = fullfile(directory, 'exports');
    
    plotEachPosition = ~true;
    
    if (reAnalyzeData)
        reAnalyze('left eye', plotEachPosition, exportsDir);
        reAnalyze('right eye', plotEachPosition, exportsDir);
        
        % Rank subjects 
        rankStrategy = 'resolution'; %  Choose from {'resolution', 'peak resolution', 'correlation coefficient 10 degs'}
        % left eye
        whichEye = 'left eye';
        rankedLeftEyeSubjectIDs = rankSubjects(exportsDir, whichEye, rankStrategy)
        % Plot ranked subject data
        plotRankedSubjects(exportsDir, rankedLeftEyeSubjectIDs, whichEye, rankStrategy);
        
        % right eye
        whichEye = 'right eye';
        rankedRightEyeSubjectIDs = rankSubjects(exportsDir, whichEye, rankStrategy)
        % Plot ranked subject data
        plotRankedSubjects(exportsDir, rankedRightEyeSubjectIDs, whichEye, rankStrategy);     
    else
        doRankAnalysis = true;
        if (doRankAnalysis)
            rankStrategy = 'resolution'; %  Choose from {'resolution', 'peak resolution', 'correlation coefficient 10 degs'}
            % Plot ranked subject data for the left eye
            whichEye = 'left eye';
            rankedLeftEyeSubjectIDs = rankSubjects(exportsDir, whichEye, rankStrategy)
            plotRankedSubjects(exportsDir, rankedLeftEyeSubjectIDs, whichEye, rankStrategy);
            % Plot ranked subject data for the right eye
            whichEye = 'right eye';
            rankedRightEyeSubjectIDs = rankSubjects(exportsDir, whichEye, rankStrategy)
            plotRankedSubjects(exportsDir, rankedRightEyeSubjectIDs, whichEye, rankStrategy);
        end
        
        plotSummary(exportsDir);
    end
end

function plotRankedSubjects(exportsDir, rankedSubjectIDs, whichEye, rankStrategy)
    dataFile = fullfile(exportsDir, sprintf('ArtalOpticsAnalysis_%s.mat', whichEye));
    load(dataFile, 'subjectPSFData', 'coneCutoffSF', 'mosaicNyquistFrequencyCPD');
    
    subjectIDs = [];
    for includedSubjectCount = 1:numel(subjectPSFData)
        d = subjectPSFData{includedSubjectCount};
        horizontalEcc = d.horizontalEcc;
        psfXCutoffSF(includedSubjectCount,:) = d.psfXCutoffSF;
        psfYCutoffSF(includedSubjectCount,:) = d.psfYCutoffSF;
        subjectIDs = cat(2, subjectIDs, d.subjectID);
    end
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
     
    subjectGroup = 0;
    for includedSubjectCount = 1:numel(rankedSubjectIDs)
        
        subjectID = rankedSubjectIDs(includedSubjectCount);
        index = find(subjectIDs == rankedSubjectIDs(includedSubjectCount));
        
        if (mod(includedSubjectCount-1,12) == 0)
            subjectGroup = subjectGroup + 1;
            hFig = figure();
            clf;
            set(hFig, 'Position', [10 10 2000 1400], 'Color', [1 1 1]);
        end
         
        r = floor((includedSubjectCount-1)/colsNum);
        r = mod(r,rowsNum)+1;
        c = mod(includedSubjectCount-1,colsNum)+1;
        subplot('Position', sv(r,c).v);
        
        plot(horizontalEcc, mosaicNyquistFrequencyCPD, 'k-', 'LineWidth', 3.0); hold on;
        plot(horizontalEcc, mosaicNyquistFrequencyCPD, 'g--', 'LineWidth', 1.0);
        if (strcmp(rankStrategy, 'correlation coefficient 10 degs'))
            plot(horizontalEcc, meanPSFXCutoffSF, 'm-', 'LineWidth', 1.5);
        end
        plot(horizontalEcc, psfXCutoffSF(index,:), 'bo-', 'MarkerFaceColor', [0.5 0.5 1], 'LineWidth', 1.5);
        plot(horizontalEcc, psfYCutoffSF(index,:), 'ro-', 'MarkerFaceColor', [1 0.5 0.5], 'LineWidth', 1.5);
        title(sprintf('%s - subjectID: %d, rank: %d/%d', whichEye, subjectID, includedSubjectCount,numel(rankedSubjectIDs)));
        
        axis 'square';
        grid on;
        set(gca, 'FontSize', 14, 'XLim', [-21 21], 'YLim', [0 70], 'XTick', -20:5:20, 'YTick', 0:5:70);
        drawnow;
        
        if (mod(includedSubjectCount-1,12) == 11)
            figTitle = sprintf('Ranked_ARTAL_Group%d_%s', subjectGroup, whichEye);
            NicePlot.exportFigToPDF(fullfile(exportsDir,sprintf('%s.pdf', figTitle)), hFig, 300);
        end
        
    end
    
end


function rankedSubjectIDs = rankSubjects(exportsDir, whichEye, rankStrategy)
    dataFile = fullfile(exportsDir, sprintf('ArtalOpticsAnalysis_%s.mat', whichEye));
    load(dataFile, 'subjectPSFData', 'coneCutoffSF', 'mosaicNyquistFrequencyCPD');
    subjectIDs = [];
    for includedSubjectCount = 1:numel(subjectPSFData)
        d = subjectPSFData{includedSubjectCount};
        horizontalEcc = d.horizontalEcc;
        psfXCutoffSF(includedSubjectCount,:) = d.psfXCutoffSF;
        psfYCutoffSF(includedSubjectCount,:) = d.psfYCutoffSF; 
        subjectIDs = cat(2, subjectIDs, d.subjectID);
    end
    
    switch (rankStrategy)
        case 'correlation coefficient 10 degs'
            % Rank according to correlation coeff in [-10 10]
            
            idx = find(abs(horizontalEcc)<=10);
            psfXCutoffSF = psfXCutoffSF(:,idx);
            psfYCutoffSF = psfYCutoffSF(:,idx);
            meanPSFXCutoffSF = mean(psfXCutoffSF,1);
            meanPSFYCutoffSF = mean(psfYCutoffSF,1);

            % Rank according to correlation coeff in [-10 10]
            r = corr(meanPSFXCutoffSF', psfXCutoffSF');
            
        case 'peak resolution'
            fovealEccX = find(horizontalEcc == 0);
            r = psfXCutoffSF(:,fovealEccX);
            
        case 'resolution'
            fovealEccX = find(horizontalEcc == 0);
            r = sqrt(psfXCutoffSF(:,fovealEccX) .* psfYCutoffSF(:,fovealEccX));
            
        otherwise
            error('Unknown rankStrategy')
    end
   
    
    [r,sortedSubjectIndices] = sort(r, 'descend');
    rankedSubjectIDs = subjectIDs(sortedSubjectIndices);
    
    hFig = figure();
    set(hFig, 'Color', [1 1 1], 'Position', [10 10 2200 800]);
    plot(1:numel(rankedSubjectIDs), r, 'bo-', 'MarkerFaceColor', [0.5 0.5 1], 'MarkerSize', 12, 'LineWidth', 1.5);
    xlabel(sprintf('subject ID (%s)', whichEye));
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
    NicePlot.exportFigToPDF(fullfile(exportsDir, sprintf('ArtalSubjectsRanked_%s', whichEye)), hFig, 300);
end


function plotSummary(exportsDir)
    dataFile = fullfile(exportsDir, sprintf('ArtalOpticsAnalysis_%s.mat', 'left eye'));
    load(dataFile, 'subjectPSFData', 'coneCutoffSF', 'mosaicNyquistFrequencyCPD');
    
    fovealEccX = [];
    for includedSubjectCount = 1:numel(subjectPSFData)
        d = subjectPSFData{includedSubjectCount};
        horizontalEcc = fliplr(d.horizontalEcc);
        if (isempty(fovealEccX))
            fovealEccX = find(horizontalEcc == 0);
        end
        psfXCutoffSF(includedSubjectCount,:) = fliplr(d.psfXCutoffSF);
        psfYCutoffSF(includedSubjectCount,:) = fliplr(d.psfYCutoffSF);
        zCoeffs(includedSubjectCount,:,:) = d.zCoeffs;
        [maxSF,idx] = max(d.psfXCutoffSF);
        bestEcc(includedSubjectCount) = horizontalEcc(idx);
        sfCutoff(includedSubjectCount) = maxSF;
    end
    


    dataFile = fullfile(exportsDir, sprintf('ArtalOpticsAnalysis_%s.mat', 'right eye'));
    load(dataFile, 'subjectPSFData', 'coneCutoffSF', 'mosaicNyquistFrequencyCPD');
    
    for includedSubjectCount = 1:numel(subjectPSFData)
        d = subjectPSFData{includedSubjectCount};
        horizontalEcc = d.horizontalEcc;
        psfXCutoffSF(includedSubjectCount+numel(subjectPSFData),:) = d.psfXCutoffSF;
        psfYCutoffSF(includedSubjectCount+numel(subjectPSFData),:) = d.psfYCutoffSF;
        [maxSF,idx] = max(d.psfXCutoffSF);
        bestEcc(includedSubjectCount+numel(subjectPSFData)) = horizontalEcc(idx);
        sfCutoff(includedSubjectCount+numel(subjectPSFData)) = maxSF;
        zCoeffs(includedSubjectCount+numel(subjectPSFData),:,:) = d.zCoeffs;
    end

    
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
       
    theZs = zCoeffs(includedSubjectCount:end,:,:);
    idx = find((horizontalEcc >= 13) & (horizontalEcc <= 18));
    theZs(:,idx,:) = nan;
  
    hFig = figure(110); clf;
    set(hFig, 'Position', [10 10 1400 975], 'Color', [1 1 1]);
    for zCoeffIndex = 4:15
        kk = zCoeffIndex-3;
        r = floor((kk-1)/colsNum);
        r = mod(r,rowsNum)+1;
        c = mod(kk-1,colsNum)+1;
        subplot('Position', sv(r,c).v);
        theZ2s = squeeze(theZs(:,:,zCoeffIndex));
        plot(horizontalEcc, theZ2s, 'b-', 'LineWidth', 1.0);
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
    
    
    % Plot distribution of Z3,Z4,Z5 at the fovea
    hFig = figure(1); clf;
    set(hFig, 'Position', [10 10 1300 510], 'Color', [1 1 1]);
    subplot(1,3,1);
    histogram(squeeze(zCoeffs(:,fovealEccX, 4)), -1:0.1:1);
    xlabel('z3 (oblique astigmatism)'); 
    axis 'square'; grid on
    set(gca, 'XTick', -1:0.2:1, 'YLim', [0 40], 'YTick', 0:10:50);
    
    subplot(1,3,2);
    histogram(squeeze(zCoeffs(:,fovealEccX, 5)), -1:0.1:1);
    xlabel('z4 (defocus)'); 
    axis 'square'; grid on
    set(gca, 'XTick', -1:0.2:1, 'YLim', [0 40], 'YTick', 0:10:50);
    
    subplot(1,3,3);
    histogram(squeeze(zCoeffs(:,fovealEccX, 6)), -1:0.1:1);
    xlabel('z5 (vertical astigmatism)'); 
    axis 'square'; grid on
    set(gca, 'XTick', -1:0.2:1, 'YLim', [0 40], 'YTick', 0:10:50);
    NicePlot.exportFigToPDF(fullfile(exportsDir,'ArtalCoeffs.pdf'), hFig, 300);
    
    
        
    % Plot the position of best resolution
    hFig = figure(2); clf;
    set(hFig, 'Position', [10 10 950 410], 'Color', [1 1 1]);
    subplot(1,2,1);
    histogram(bestEcc, (-20:1:20)-0.5);
    text(14, 3, 'OD', 'FontSize', 12, 'FontWeight', 'Bold');
    xlabel('horizontal eccentricity of peak SF cutoff (degs)');
    ylabel('count');  
    set(gca, 'FontSize', 12, 'XTick', -20:5:20, 'XLim', [-21 21]);
    axis 'square'
    
    subplot(1,2,2);
    plot(bestEcc, sfCutoff, 'k.');
    xlabel('horizontal eccentricity of peak SF cutoff(degs)');
    ylabel('peak SF cutoff (-3dB) (c/deg)');
    set(gca, 'FontSize', 12, 'XTick', -20:5:20, 'XLim', [-21 21], 'YLim', [0 65]);
    text(14, 2, 'OD', 'FontSize', 12, 'FontWeight', 'Bold');
    axis 'square'
    NicePlot.exportFigToPDF(fullfile(exportsDir,'ArtalBestEcc.pdf'), hFig, 300);

    
    
    hFig = figure(1);
    for vIndex = 1:size(mosaicNyquistFrequencyCPD,1)
        hFig = figure(1);
        clf;
        set(hFig, 'Color', [1 1 1], 'Position', [10 10 700 700]);
        h1 = plot(horizontalEcc, mosaicNyquistFrequencyCPD(vIndex,:), 'k-', 'LineWidth', 3);
        hold on
        plot(horizontalEcc, mosaicNyquistFrequencyCPD(vIndex,:), 'g--', 'LineWidth', 1.5);
        h3 = plot(horizontalEcc, mean(psfXCutoffSF,1), 'bo-', 'LineWidth', 1.5, 'MarkerSize', 10, 'MarkerFaceColor', [0.5 0.5 1]);
        h4 = plot(horizontalEcc, mean(psfYCutoffSF,1), 'ro-', 'LineWidth', 1.5, 'MarkerSize', 10, 'MarkerFaceColor', [1 0.5 0.5]);
        
       % if (verticalEcc(vIndex) == 0)
        legend([h1 h3 h4], {'mosaic Nyquist freq.', 'PSF_x (mean)', 'PSF_y (mean)'});
        xlabel('horizontal eccentricity (deg)');
        %end
        
        %ylabel('spatial frequency cutoff, -15dB (c/deg)');
        %title(sprintf('%2.0f deg', verticalEcc(vIndex)));
        axis 'square';  grid 'on'; box 'off'
        set(gca, 'XLim', [-20 20], 'YLim', [0 70], 'XTick', -20:5:20, 'YTick', 0:5:100,  'FontSize', 16);
        set(gca, 'Color', 'none', 'XColor', 'none', 'YColor', 'none', 'LineWidth', 1.0);
       % if (verticalEcc(vIndex) == 0)
            set(gca, 'Color', 'none', 'XColor', [0.3 0.3 0.3], 'YColor', [0.3 0.3 0.3]);
        %end
        
        NicePlot.exportFigToPDF(fullfile(exportsDir,sprintf('ArtalMeanAcrossHorizontalEccentricity.pdf')), hFig, 300);
    end
    
    
    hFig = figure(3);
    clf;
    set(hFig, 'Color', [1 1 1], 'Position', [10 10 1500 710]);
    ax = subplot(1,2,1);
    h1 = plot(ax,horizontalEcc, mosaicNyquistFrequencyCPD, 'k-', 'LineWidth', 3);
    hold(ax, 'on');
    h2 = plot(ax,horizontalEcc, mosaicNyquistFrequencyCPD, 'g--', 'LineWidth', 1.5);
    h3 = plot(ax,horizontalEcc, mean(psfXCutoffSF,1), 'b-', 'LineWidth', 1.5);
    %h4 = plot(ax,horizontalEcc, prctile(psfXCutoffSF,95,1), 'r--', 'LineWidth', 1.5);
    %h5 = plot(ax,horizontalEcc, prctile(psfXCutoffSF,5,1), 'r--', 'LineWidth', 1.5);
    h5 = shadeAreaBetweenCyrves(ax, horizontalEcc, prctile(psfXCutoffSF,95,1), prctile(psfXCutoffSF,5,1), [0.5 0.5 1], 0.5);
    legend([h1 h3 h5], {'mosaic Nyquist freq.', 'PSF cutoff (mean)', 'PSF cutoff (5%, 95%)'});
    xlabel('horizontal eccentricity (deg)');
    ylabel('spatial frequency cutoff, -15dB (c/deg)');
    title('high-resolution axis');
    axis 'square';  grid 'on'; box 'off'
    set(gca, 'XLim', [-20 20], 'YLim', [0 70], 'XTick', -20:5:20, 'YTick', 0:5:100,  'FontSize', 16);
    set(gca, 'XColor', [0.3 0.3 0.3], 'YColor', [0.3 0.3 0.3], 'LineWidth', 1.0)
    
    ax = subplot(1,2,2);
    h1 = plot(ax,horizontalEcc, mosaicNyquistFrequencyCPD, 'k-', 'LineWidth', 3);
    hold(ax, 'on');
    h2 = plot(ax,horizontalEcc, mosaicNyquistFrequencyCPD, 'g--', 'LineWidth', 1.5);
    h3 = plot(ax,horizontalEcc, mean(psfYCutoffSF,1), 'r-', 'LineWidth', 1.5);
    %h4 = plot(horizontalEcc, prctile(psfYCutoffSF,95,1), 'b--', 'LineWidth', 1.0);
    %h5 = plot(horizontalEcc, prctile(psfYCutoffSF,5,1), 'b--', 'LineWidth', 1.0);
    h5 = shadeAreaBetweenCyrves(ax, horizontalEcc, prctile(psfYCutoffSF,95,1), prctile(psfYCutoffSF,5,1), [1 0.5 0.5], 0.5);
    legend([h1 h3 h5], {'mosaic Nyquist freq.', 'PSF cutoff (mean)', 'PSF cutoff (5%, 95%)'});
    xlabel('horizontal eccentricity (deg)');
    ylabel('spatial frequency cutoff, -15dB (c/deg)');
    title('low-resolution axis');
    axis 'square';  grid 'on'; box 'off'
    set(gca, 'XLim', [-20 20], 'YLim', [0 70], 'XTick', -20:5:20, 'YTick', 0:5:100, 'FontSize', 16);
    set(gca, 'XColor', [0.3 0.3 0.3], 'YColor', [0.3 0.3 0.3], 'LineWidth', 1.0);
    NicePlot.exportFigToPDF(fullfile(exportsDir,sprintf('ArtalAcrossHorizontalEccentricity.pdf')), hFig, 300);
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

function reAnalyze(whichEye, plotEachPosition, exportsDir)
    opticsParams = struct(...
        'zernikeDataBase',  'Artal2012', ...
        'subjectID',1, ...
        'subtractCentralRefraction', false, ...
        'zeroCenterPSF', true, ...
        'flipPSFUpsideDown', true, ...
        'whichEye', 'right eye', ...
        'pupilDiameterMM', 3);
   

    horizontalEcc = -20:20;
    whichEyes = {whichEye};
    
    
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
     
    subjectGroup = 0;
    includedSubjects = ArtalOptics.subjectsWithFullDataSets(whichEyes{1});
    
    eyeIndex = 1;
    
    for includedSubjectCount = 1:numel(includedSubjects)
        
        subjectID = includedSubjects(includedSubjectCount);
        if (mod(includedSubjectCount-1,12) == 0)
            subjectGroup = subjectGroup + 1;
            hFig = figure(subjectGroup);
            clf;
            set(hFig, 'Position', [10 10 2000 1400]);
        end
        
        
        opticsParams.subjectID = subjectID;
        opticsParams.whichEye = whichEyes{eyeIndex};
        opticsParams.subtractCentralRefraction = ArtalOptics.constants.subjectRequiresCentralRefractionCorrection(opticsParams.whichEye, opticsParams.subjectID);
        
        if (plotEachPosition)
            videoFileName = fullfile(exportsDir,sprintf('%s_Subject%d_%s_Pupil%3.2fmm_SubtractCentralRefraction', ...
                opticsParams.zernikeDataBase,  opticsParams.subjectID, opticsParams.whichEye, opticsParams.pupilDiameterMM));
            videoOBJ = VideoWriter(videoFileName, 'MPEG-4');
            videoOBJ.FrameRate = 30;
            videoOBJ.Quality = 100;
            videoOBJ.open();
        end
        
        for iEcc = 1:numel(horizontalEcc)
            
            % Update params
            mosaicEcc = [horizontalEcc(iEcc) 0];
            
            
            % Compute PSFimage and cone image (single aperture)
            [psfImage, coneImage, NyquistFrequency, psfSupportArcMin, theZCoeffs] = psfAndConeImage(mosaicEcc, opticsParams);
            if (isempty(psfImage))
                fprintf('No data for ecc (x:%2.0f)\n', horizontalEcc(iEcc));
                continue;
            end
            mosaicNyquistFrequencyCPD(iEcc) = NyquistFrequency;
            zCoeffs(iEcc,:) = theZCoeffs; 
            
            
            
            % Compute power spectra
            [psfImageSpectrum, psfImageSpectrumRoatated, coneMosaicImageSpectrum, ...
                spectralSupportCyclesPerDegree, spectralSupportCyclesPerDegreePositive, coneSpectrumSlice, ...
                spectralSupportCyclesPerDegreePositiveX, psfSpectrumSliceX, ...
                spectralSupportCyclesPerDegreePositiveY, psfSpectrumSliceY, ...
                effectivePSFSpectrumSliceX, effectivePSFSpectrumSliceY, ...
                coneCutoffSF(iEcc), psfXCutoffSF(iEcc), psfYCutoffSF(iEcc)] = ...
                computePSFAndConePowerSpectra(psfImage, coneImage, psfSupportArcMin);
            
            
            if (plotEachPosition)
                hFigVideo = visualizeSpectralAnalysis(1000,horizontalEcc(iEcc), 0, ...
                    psfSupportArcMin, psfImage, coneImage, ...
                    spectralSupportCyclesPerDegree, spectralSupportCyclesPerDegreePositive, ...
                    psfImageSpectrum, psfImageSpectrumRoatated,  ...
                    coneMosaicImageSpectrum, coneSpectrumSlice, ...
                    spectralSupportCyclesPerDegreePositiveX, psfSpectrumSliceX, effectivePSFSpectrumSliceX, ...
                    spectralSupportCyclesPerDegreePositiveY, psfSpectrumSliceY, effectivePSFSpectrumSliceY, ...
                    coneCutoffSF(iEcc), psfXCutoffSF(iEcc), psfYCutoffSF(iEcc), ...
                    mosaicNyquistFrequencyCPD(iEcc));
                % Add frame to video
                drawnow;
                videoOBJ.writeVideo(getframe(hFigVideo));
            end
     
        end
        
        % Data for exporting
        subjectPSFData{includedSubjectCount} = struct(...
            'horizontalEcc', horizontalEcc, ...
            'verticalEcc', 0, ...
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
        
        r = floor((includedSubjectCount-1)/colsNum);
        r = mod(r,rowsNum)+1;
        c = mod(includedSubjectCount-1,colsNum)+1;
        fprintf(' subject %d (index:%d (r:%d, c:%d)), subtract central refraction: %d\n', subjectID, includedSubjectCount, r,c,opticsParams.subtractCentralRefraction);
        
        hFig = figure(subjectGroup);
        subplot('Position', sv(r,c).v);
        
        h1 = plot(horizontalEcc, mosaicNyquistFrequencyCPD, 'k-', 'LineWidth', 3);
        hold on;
        plot(horizontalEcc, mosaicNyquistFrequencyCPD, 'g--', 'LineWidth', 1);
        h2 = plot(horizontalEcc, psfXCutoffSF, 'bo-', 'LineWidth', 1.5, 'MarkerSize', 8, 'MarkerFaceColor', [0.5 0.5 1]);
        h3 = plot(horizontalEcc, psfYCutoffSF, 'ro-', 'LineWidth', 1.5, 'MarkerSize', 8, 'MarkerFaceColor', [1 0.5 0.5]);
        
        if (r == 3)
            xlabel('horizontal eccentricity (deg)');
        end
        if (c == 1)
            ylabel('spatial frequency cutoff (-3dB), c/deg');
        end
        if (opticsParams.subtractCentralRefraction)
            title(sprintf('subject ID: %d (-central refraction)', opticsParams.subjectID));
        else
            title(sprintf('subject ID: %d (original)', opticsParams.subjectID));
        end
        legend([h1 h2 h3], {'mosaic Nyquist freq.', 'psf_x cutoff (-15dB)', 'psf_y cutoff (-15dB)'});
        
        axis 'square';
        grid on;
        set(gca, 'FontSize', 14, 'XLim', [-21 21], 'YLim', [0 70], 'XTick', -20:5:20, 'YTick', 0:5:70);
        drawnow;
        
        if (mod(includedSubjectCount-1,12) == 11)
            figTitle = sprintf('%s_Group%d_%s_Pupil%3.2fmm', ...
                    opticsParams.zernikeDataBase, subjectGroup, opticsParams.whichEye, opticsParams.pupilDiameterMM);
            NicePlot.exportFigToPDF(fullfile(exportsDir,sprintf('%s.pdf', figTitle)), hFig, 300);
        end
        

    end
    
    % Export data
    dataFile = sprintf(fullfile(exportsDir,sprintf('ArtalOpticsAnalysis_%s.mat', whichEye)));
    save(dataFile, 'subjectPSFData', 'coneCutoffSF', 'mosaicNyquistFrequencyCPD');
    
end
