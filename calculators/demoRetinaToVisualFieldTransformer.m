function demoRetinaToVisualFieldTransformer()

    
    xFormer = RetinaToVisualFieldTransformer('ZernikeDataBase', 'Artal2012');

    
    analyzedRetinalQuadrant = RetinaToVisualFieldTransformer.nasalRetinaQuadrant;
    subjectRankingEye = RetinaToVisualFieldTransformer.rightEye;
    analyzedEye = RetinaToVisualFieldTransformer.rightEye;
    pupilDiameterMM = 3.0;

    % All subjects
    examinedSubjectRankOrders = 1:RetinaToVisualFieldTransformer.ArtalSubjectsNum;

    % Only the first 30
    examinedSubjectRankOrders = 1:38;
    % Remove some subjects which increase the
    examinedSubjectRankOrders = setdiff(examinedSubjectRankOrders, [5 16 20 23 31 34 36 37]);


    % Get the range of eccentricities for this retinal quadrant
    maxEccDegs = 30;
    [horizontalEccDegs, verticalEccDegs, eccDegsForPlotting] = ...
         RetinaToVisualFieldTransformer.eccentricitiesForQuadrant(...
                analyzedRetinalQuadrant, analyzedEye, maxEccDegs);

    % Figures and video
    generateFigsAndVideos = false;

    % Allocate memory
    visualConeCharacteristicRadiusDegs = nan(numel(examinedSubjectRankOrders), numel(horizontalEccDegs));
    anatomicalConeCharacteristicRadiusDegs = nan(1, numel(horizontalEccDegs));
    conesNumInPatch = nan(1, numel(horizontalEccDegs));

    % Analyze all examined subjects
    for iSubj = 1:numel(examinedSubjectRankOrders)

        % Get the subject ID
        subjectRankOrder = examinedSubjectRankOrders(iSubj);
        subjID = xFormer.subjectWithRankInEye(subjectRankOrder, subjectRankingEye);

        % Assemble filename
        dataFileName = sprintf('%s_SubjectID%d_%s_%s_PupilDiam%2.2fMM', ...
            xFormer.ZernikeDataBase, subjID, analyzedEye, upper(strrep(analyzedRetinalQuadrant, ' ', '_')), pupilDiameterMM);

        if (generateFigsAndVideos)
            p = getpref('ISETMacaque');
            videoFName = fullfile(p.generatedDataDir, 'coneApertureBackingOut', dataFileName);
    
            videoOBJ = VideoWriter(videoFName, 'MPEG-4');
            videoOBJ.FrameRate = 30;
            videoOBJ.Quality = 100;
            videoOBJ.open();
        else
            videoOBJ = [];
        end

        
        parfor iEcc = 1:numel(horizontalEccDegs)
            if (generateFigsAndVideos)
                hFig = figure(10); clf;
                set(hFig, 'Color', [1 1 1], 'Position', [10 10 1600 400]);
            else
                hFig = [];
            end
            
            % Analyze effect of optics at this eccentricity
            eccDegs = [horizontalEccDegs(iEcc) verticalEccDegs(iEcc)];
            
            dStruct = xFormer.estimateConeCharacteristicRadiusInVisualSpace(...
                analyzedEye, eccDegs, subjID, pupilDiameterMM, dataFileName, ...
                'anatomicalConeCharacteristicRadiusDegs', [], ...
                'hFig', hFig, ...
                'videoOBJ', videoOBJ);

            % Save results
            conesNumInPatch(iEcc) = dStruct.conesNumInRetinalPatch;
            anatomicalConeCharacteristicRadiusDegs(iEcc) = dStruct.anatomicalConeCharacteristicRadiusDegs;
            visualConeCharacteristicRadiusDegs(iSubj,iEcc) = dStruct.visualConeCharacteristicRadiusDegs;
        end %iEcc

        if (~isempty(videoOBJ))
            videoOBJ.close();
        end


        % Summary data

        % Current subject data
        visualToAnatomicalRcRatio = visualConeCharacteristicRadiusDegs(iSubj,:) ./ anatomicalConeCharacteristicRadiusDegs(1,:);

        p = getpref('ISETMacaque');
        pdfFileName = fullfile(p.generatedDataDir, 'coneApertureBackingOut', sprintf('%s_Summary.pdf',dataFileName));

        plotData(100+subjectRankOrder, pdfFileName, eccDegsForPlotting, conesNumInPatch,  ...
            visualToAnatomicalRcRatio, [], ...
            analyzedRetinalQuadrant, analyzedEye);

        % Mean over subjects data
        visualToAnatomicalRcRatioMeanOverSubjects = mean(visualConeCharacteristicRadiusDegs,1, 'omitnan') ./ anatomicalConeCharacteristicRadiusDegs(1,:);
        visualToAnatomicalRcRatioSTDOverSubjects = std(visualConeCharacteristicRadiusDegs,0,1, 'omitnan') ./ anatomicalConeCharacteristicRadiusDegs(1,:);
        
        
        pdfFileName = fullfile(p.generatedDataDir, 'coneApertureBackingOut', sprintf('%s_Summary_%dSubjects.pdf',dataFileName, iSubj));
        plotData(1000+subjectRankOrder, pdfFileName, eccDegsForPlotting, conesNumInPatch, ...
            visualToAnatomicalRcRatioMeanOverSubjects, visualToAnatomicalRcRatioSTDOverSubjects, ...
            analyzedRetinalQuadrant, analyzedEye);

    end %iSubj
end

function plotData(figNo, pdfFilename, eccDegsForPlotting, conesNumInPatch,  visualToAnatomicalRcRatio, visualToAnatomicalRcRatioSTD, analyzedRetinalQuadrant, analyzedEye)
    
    % Plot data
    hFig = figure(figNo); clf;
    set(hFig, 'Color', [1 1 1], 'Position', [100 100 1250 520]);
    ax = subplot(1,2,1);
    plot(ax, eccDegsForPlotting, conesNumInPatch, 'rs-', 'MarkerSize', 10, 'MarkerFaceColor', [1 0.5 0.5], 'LineWidth', 1.5); 
    set(ax, 'XLim', [-1 31], 'XTick', 0:2:30);
    if (strcmp(analyzedRetinalQuadrant, RetinaToVisualFieldTransformer.temporalRetinaQuadrant))
        set(ax, 'XDir', 'reverse')
    end
    xtickangle(ax, 0);
    grid(ax, 'on');
    legend(ax,{analyzedEye}, 'Location', 'northeast')
    ylabel(ax,'# of cones/deg2')
    xlabel(ax, sprintf('eccentricity, %s (degs)', analyzedRetinalQuadrant));
    set(ax, 'FontSize', 16)

    ax = subplot(1,2,2);
    if (~isempty(visualToAnatomicalRcRatioSTD))
        if (strcmp(analyzedRetinalQuadrant, RetinaToVisualFieldTransformer.nasalRetinaQuadrant))
            n1 = 14;
            patch1 = patchCoordsFromXYmeanYstd(eccDegsForPlotting(1:n1), visualToAnatomicalRcRatio(1:n1), visualToAnatomicalRcRatioSTD(1:n1));

            n2 = 19;
            patch2 = patchCoordsFromXYmeanYstd(eccDegsForPlotting(n2:end), visualToAnatomicalRcRatio(n2:end), visualToAnatomicalRcRatioSTD(n2:end));
        else
            patch1 = patchCoordsFromXYmeanYstd(eccDegsForPlotting, visualToAnatomicalRcRatio, visualToAnatomicalRcRatioSTD);
            patch2 = [];
        end

        patch(ax, patch1.x, patch1.y, 0*patch1.y, [1 0 0], 'FaceColor', [1 0.5 0.5], 'FaceAlpha', 0.5);
        if (~isempty(patch2))
            patch(ax, patch2.x, patch2.y, 0*patch2.y, [1 0 0], 'FaceColor', [1 0.5 0.5], 'FaceAlpha', 0.5);
        end
        hold(ax, 'on');
    end

    plot(ax, eccDegsForPlotting, visualToAnatomicalRcRatio, 'ro-', 'MarkerSize', 10, 'LineWidth', 1.5, 'MarkerFaceColor', [1 0.5 0.5]); hold on
    
    set(ax, 'XLim', [-1 31], 'XTick', 0:2:30, 'YLim', [0 16], 'YTick', 0:2:20);
    xtickangle(ax, 0);
    if (strcmp(analyzedRetinalQuadrant, RetinaToVisualFieldTransformer.temporalRetinaQuadrant))
        set(ax, 'XDir', 'reverse')
    end

    grid(ax, 'on');
    legend(ax,{analyzedEye}, 'Location', 'northeast')
    ylabel(ax,'visual cone Rc/anatomical cone Rc')
    xlabel(ax, sprintf('eccentricity, %s (degs)', analyzedRetinalQuadrant));
    set(ax, 'FontSize', 16)
    drawnow;

    NicePlot.exportFigToPDF(pdfFilename, hFig, 300);

end

function thePatch = patchCoordsFromXYmeanYstd(x, yMean, ySTD)
    thePatch.x(1) = x(1);
    thePatch.y(1) = yMean(1)-ySTD(1);

    thePatch.x(1+(1:numel(x))) = x;
    thePatch.y(1+(1:numel(x))) = yMean + ySTD;

    thePatch.x(numel(thePatch.x)+1) = thePatch.x(numel(thePatch.x));
    thePatch.y(numel(thePatch.y)+1) = yMean(end) - ySTD(end);

    thePatch.x(numel(thePatch.x)+(1:numel(x))) = fliplr(x);
    thePatch.y(numel(thePatch.y)+(1:numel(x))) = fliplr(yMean - ySTD);
end
