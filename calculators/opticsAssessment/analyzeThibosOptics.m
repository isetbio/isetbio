function analyzeThibosOptics(reAnalyzeData)

    % Get directory
    [directory,~] = fileparts(which(mfilename()));
    exportsDir = fullfile(directory, 'exports');

    if (reAnalyzeData)
        reAnalyze('left eye', exportsDir);
        reAnalyze('right eye', exportsDir);
    else
        doRankAnalysis = true;
        if (doRankAnalysis)
            rankStrategy = 'resolution'; %  Choose from {'resolution', 'peak resolution', 'correlation coefficient 10 degs'}
           
            % Plot ranked subject data for the left eye
            whichEye = 'left eye';
            rankedLeftEyeSubjectIDs = rankSubjects(exportsDir, whichEye, rankStrategy)
            pause
            plotSummary(exportsDir,whichEye);

            % Plot ranked subject data for the left eye
            whichEye = 'right eye';
            rankedRightEyeSubjectIDs = rankSubjects(exportsDir, whichEye, rankStrategy)
            plotSummary(exportsDir,whichEye);
        end

        
    end
end

function rankedSubjectIDs = rankSubjects(exportsDir, whichEye, rankStrategy)

    dataFile = fullfile(exportsDir, sprintf('ThibosOpticsAnalysis_%s.mat', whichEye));
    load(dataFile, 'subjectPSFData', 'coneCutoffSF', 'mosaicNyquistFrequencyCPD');
    
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
        case 'resolution'
            ylabel('foveal resolution (c/deg)');
            set(gca, 'YLim', [0 65]);
    end
    set(gca, 'FontSize', 16, 'XLim',[0 numel(rankedSubjectIDs)+1], 'XTick', 1:numel(rankedSubjectIDs), 'XTickLabel', rankedSubjectIDs); 
    grid on
    NicePlot.exportFigToPDF(fullfile(exportsDir, sprintf('ThibosSubjectsRanked_%s', whichEye)), hFig, 300);


end

function plotSummary(exportsDir,whichEye)
    dataFile = sprintf(fullfile(exportsDir,sprintf('ThibosOpticsAnalysis_%s.mat', whichEye)));  
    load(dataFile, 'subjectPSFData', 'coneCutoffSF', 'mosaicNyquistFrequencyCPD');
    for includedSubjectCount = 1:numel(subjectPSFData)
        d = subjectPSFData{includedSubjectCount};
        psfXCutoffSF(includedSubjectCount) = d.psfXCutoffSF;
        psfYCutoffSF(includedSubjectCount) = d.psfYCutoffSF;
    end
    
    pause
    
    hFig = figure();
    clf;
    set(hFig, 'Color', [1 1 1], 'Position', [10 10 740 710]);
    scatter(psfXCutoffSF, psfYCutoffSF, 14*14, 'o', ...
        'MarkerFaceColor', [1 0.5 0.5], 'MarkerEdgeColor', [0.3 0.3 0.3], 'MarkerFaceAlpha', 0.5, 'LineWidth', 1.5);
    hold on;
    plot([0 100], mosaicNyquistFrequencyCPD*[1 1], 'k-', 'LineWidth', 3.0);
    plot([0 100], mosaicNyquistFrequencyCPD*[1 1], 'g-', 'LineWidth', 1.0);
    plot(mosaicNyquistFrequencyCPD*[1 1], [0 100], 'k-', 'LineWidth', 3.0);
    plot(mosaicNyquistFrequencyCPD*[1 1], [0 100], 'g--', 'LineWidth', 1.0);
    legend({'Thibos 2009'});
    title(sprintf('Thibos (%s)', whichEye));
    xlabel('spatial frequency (high-resolution axis)');
    ylabel('spatial frequency (low-resolution axis)');
    axis 'square';  grid 'on'
    set(gca, 'XLim', [0 100], 'YLim', [0 100], 'FontSize', 16);
    set(gca, 'XColor', [0.3 0.3 0.3], 'YColor', [0.3 0.3 0.3], 'LineWidth', 1.0)
end


function reAnalyze(whichEye, exportsDir)
    opticsParams = struct(...
        'zernikeDataBase',  'Thibos2002', ...
        'subjectID',1, ...
        'subtractCentralRefraction', true, ...
        'zeroCenterPSF', true, ...
        'flipPSFUpsideDown', true, ...
        'whichEye', 'right eye', ...
        'pupilDiameterMM', 3);
    
    includedSubjects = 1:70;
    whichEyes = {whichEye};
    
    plotEachPosition = true;
    subtractCentralRefraction = false;
    
    eyeIndex = 1;
    
    mosaicEcc = [0 0];
    
    if (plotEachPosition)
            videoFileName = sprintf('%s_%s_Pupil%3.2fmm_SubtractCentralRefraction', ...
                opticsParams.zernikeDataBase,  opticsParams.whichEye, opticsParams.pupilDiameterMM);
            videoOBJ = VideoWriter(videoFileName, 'MPEG-4');
            videoOBJ.FrameRate = 30;
            videoOBJ.Quality = 100;
            videoOBJ.open();
    end
        
    subjectGroup = 0;
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
        opticsParams.subtractCentralRefraction = (subtractCentralRefraction==1);
        
        % Compute PSFimage and cone image (single aperture)
        [psfImage, coneImage, NyquistFrequency, psfSupportArcMin, zCoeffs] = psfAndConeImage(mosaicEcc, opticsParams);
            
        if (isempty(psfImage))
            fprintf('No data for ecc (x:%2.0f)\n', horizontalEcc(iEcc));
            continue;
        end
        mosaicNyquistFrequencyCPD = NyquistFrequency;

        % Compute power spectra
        [psfImageSpectrum, psfImageSpectrumRoatated, coneMosaicImageSpectrum, ...
            spectralSupportCyclesPerDegree, spectralSupportCyclesPerDegreePositive, coneSpectrumSlice, ...
            spectralSupportCyclesPerDegreePositiveX, psfSpectrumSliceX, ...
            spectralSupportCyclesPerDegreePositiveY, psfSpectrumSliceY, ...
            effectivePSFSpectrumSliceX, effectivePSFSpectrumSliceY, ...
            coneCutoffSF(includedSubjectCount), psfXCutoffSF(includedSubjectCount), psfYCutoffSF(includedSubjectCount)] = ...
            computePSFAndConePowerSpectra(psfImage, coneImage, psfSupportArcMin);

        % Data for exporting
        subjectPSFData{includedSubjectCount} = struct(...
            'horizontalEcc', 0, ...
            'verticalEcc', 0, ...
            'subjectID', opticsParams.subjectID, ...
            'zCoeffs', zCoeffs, ...
            'whichEye', opticsParams.whichEye, ...
            'pupilDiamMM', opticsParams.pupilDiameterMM, ...
            'psfXCutoffSF', psfXCutoffSF(includedSubjectCount), ...
            'psfYCutoffSF', psfYCutoffSF(includedSubjectCount) ...
            );
        
        if (plotEachPosition)
            hFigVideo = visualizeSpectralAnalysis(1000, 0, 0, ...
                psfSupportArcMin, psfImage, coneImage, ...
                spectralSupportCyclesPerDegree, spectralSupportCyclesPerDegreePositive, ...
                psfImageSpectrum, psfImageSpectrumRoatated,  ...
                coneMosaicImageSpectrum, coneSpectrumSlice, ...
                spectralSupportCyclesPerDegreePositiveX, psfSpectrumSliceX, effectivePSFSpectrumSliceX, ...
                spectralSupportCyclesPerDegreePositiveY, psfSpectrumSliceY, effectivePSFSpectrumSliceY, ...
                coneCutoffSF(includedSubjectCount), psfXCutoffSF(includedSubjectCount), psfYCutoffSF(includedSubjectCount), ...
                mosaicNyquistFrequencyCPD);
            % Add frame to video
            drawnow;
            videoOBJ.writeVideo(getframe(hFigVideo));
        end
    end
    
    if (plotEachPosition)
        videoOBJ.close();
    end
        
    % Export data
    dataFile = sprintf(fullfile(exportsDir,sprintf('ThibosOpticsAnalysis_%s.mat', whichEye)));
    save(dataFile, 'subjectPSFData', 'coneCutoffSF', 'mosaicNyquistFrequencyCPD');
end


            
        