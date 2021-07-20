function contrastPolansToArtal

     % Get directory
    [directory,~] = fileparts(which(mfilename()));
    exportsDir = fullfile(directory, 'exports');
    
    [PolansZcoeffs, horizontalEccPolans] = loadPolansZCoeffs(exportsDir);
    [ArtalZcoeffs, horizontalEccArtal] = loadArtalZCoeffs(exportsDir);
    
    rowsNum = 3;
    colsNum = 4;
    sv = NicePlot.getSubPlotPosVectors(...
           'colsNum', colsNum, ...
           'rowsNum', rowsNum, ...
           'heightMargin',  0.06, ...
           'widthMargin',    0.02, ...
           'leftMargin',     0.02, ...
           'rightMargin',    0.00, ...
           'bottomMargin',   0.05, ...
           'topMargin',      0.01); 
       
    hFig = figure(111); clf;
    set(hFig, 'Position', [10 10 1400 975], 'Color', [1 1 1]);
    for zCoeffIndex = 4:15
        kk = zCoeffIndex-3;
        r = floor((kk-1)/colsNum);
        r = mod(r,rowsNum)+1;
        c = mod(kk-1,colsNum)+1;
        subplot('Position', sv(r,c).v);
        
        Z1 = squeeze(ArtalZcoeffs(:,:,zCoeffIndex));
        Z2 = squeeze(PolansZcoeffs(:,:,zCoeffIndex));
        maxZ = max([max(abs(Z1(:))) max(abs(Z2(:)))]);
        
        plot(horizontalEccArtal, Z1, 'ro-', 'MarkerFaceColor', [1 0.5 0.5], 'LineWidth', 1.0);
        hold on;
        plot(horizontalEccPolans, Z2, 'b-', 'MarkerFaceColor', [0.5 0.5 1], 'LineWidth', 1.0);
        ylabel(sprintf('Z%d', zCoeffIndex-1));
        axis 'square';
        grid on;
        maxZamp = max([0.35 maxZ]);
        set(gca, 'FontSize', 14, 'YLim', maxZamp*[-1 1], 'XLim', [-21 21], 'XTick', -20:5:20);
        if (zCoeffIndex == 5)
            ylabel(sprintf('Z%d (defocus)', zCoeffIndex-1));
        elseif (zCoeffIndex == 4)
            ylabel(sprintf('Z%d (oblique astigmatism)', zCoeffIndex-1));
        elseif (zCoeffIndex == 6)
            ylabel(sprintf('Z%d (vertical astigmatism)', zCoeffIndex-1));
        end
        xlabel('eccentricity (degs)');
    end
    
    NicePlot.exportFigToPDF(fullfile(exportsDir,'PolansVsArtalZCoeffsHorizontalMeridian.pdf'), hFig, 300);
    
end

function [ArtalZcoeffs, horizontalEccArtal] = loadArtalZCoeffs(exportsDir)

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

    theZs = zCoeffs(includedSubjectCount:end,:,:);
    idx = find((horizontalEcc >= 13) & (horizontalEcc <= 18));
    theZs(:,idx,:) = nan;
    
    ArtalZcoeffs = theZs;
    horizontalEccArtal = horizontalEcc;
end

function [PolansZcoeffs, horizontalEccPolans] = loadPolansZCoeffs(exportsDir)
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
    
   
       
    theZs = squeeze(zCoeffs(:,fovealEccY,:,:));
    idx = find((horizontalEcc >= 13) & (horizontalEcc <= 18));
    theZs(:,idx,:) = nan;
  
    PolansZcoeffs = theZs;
    horizontalEccPolans = horizontalEcc;
end
