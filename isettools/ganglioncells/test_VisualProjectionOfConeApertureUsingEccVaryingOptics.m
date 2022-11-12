function test_VisualProjectionOfConeApertureUsingEccVaryingOptics()

    % Intantiate RetinaToVisualFieldTrasformer with the Artal database
    xFormer = RetinaToVisualFieldTransformer('ZernikeDataBase', 'Artal2012');

    % Analyze along the temporal meridian
    analyzedRetinaMeridian = 'temporal meridian';
    subjectRankingEye = 'right eye';
    analyzedEye = 'right eye';
    pupilDiameterMM = 3.0;

    % All subjects
    examinedSubjectRankOrders = 1:RetinaToVisualFieldTransformer.ArtalSubjectsNum;

    % Only the first 30
    examinedSubjectRankOrders = 1:38;
    % Remove some subjects which increase the
    examinedSubjectRankOrders = setdiff(examinedSubjectRankOrders, [5 16 20 23 31 34 36 37]);


    % Get the horizontal and vertical eccs corresponding to the target radial ecc at the
    % analyzed retinal meridian and eye
    radialEccDegs = 0:1:30;
    [horizontalEccDegs, verticalEccDegs] = cMosaic.eccentricitiesForRetinaMeridianInEye(...
            radialEccDegs, analyzedRetinaMeridian, analyzedEye);

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
            xFormer.ZernikeDataBase, subjID, analyzedEye, upper(strrep(analyzedRetinaMeridian, ' ', '_')), pupilDiameterMM);

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
            visualToAnatomicalRadiiRatios(iSubj,iEcc) = visualConeCharacteristicRadiusDegs(iSubj,iEcc) / anatomicalConeCharacteristicRadiusDegs(iEcc);
        end %iEcc

        if (~isempty(videoOBJ))
            videoOBJ.close();
        end


        % Summary data

        % Current subject data
        visualToAnatomicalRcRatio = visualConeCharacteristicRadiusDegs(iSubj,:) ./ anatomicalConeCharacteristicRadiusDegs(1,:);

        p = getpref('ISETMacaque');
        pdfFileName = fullfile(p.generatedDataDir, 'coneApertureBackingOut', sprintf('%s_Summary.pdf',dataFileName));

        plotData(100+subjectRankOrder, pdfFileName, radialEccDegs, conesNumInPatch,  ...
            visualToAnatomicalRcRatio, [], ...
            analyzedRetinaMeridian, analyzedEye);

        % Mean over subjects data
        

        visualToAnatomicalRcRatioMeanOverSubjects = mean(visualToAnatomicalRadiiRatios,1, 'omitnan'); % mean(visualConeCharacteristicRadiusDegs,1, 'omitnan') ./ anatomicalConeCharacteristicRadiusDegs(1,:);
        visualToAnatomicalRcRatioPrcTilesOverSubjects = prctile(visualToAnatomicalRadiiRatios, [5 95 25 75], 1); %; %std(visualConeCharacteristicRadiusDegs,0,1, 'omitnan') ./ anatomicalConeCharacteristicRadiusDegs(1,:);
        

        pdfFileName = fullfile(p.generatedDataDir, 'coneApertureBackingOut', sprintf('%s_Summary_%dSubjects.pdf',dataFileName, iSubj));
        plotData(1000+subjectRankOrder, pdfFileName, radialEccDegs, conesNumInPatch, ...
            visualToAnatomicalRcRatioMeanOverSubjects, visualToAnatomicalRcRatioPrcTilesOverSubjects, ...
            analyzedRetinaMeridian, analyzedEye);

    end %iSubj
end

function plotData(figNo, pdfFilename, radialEccDegs, conesNumInPatch,  visualToAnatomicalRcRatio, visualToAnatomicalRcRatioPrcTilesOverSubjects, analyzedRetinaMeridian, analyzedEye)
    
    % Plot data
    hFig = figure(figNo); clf;
    set(hFig, 'Color', [1 1 1], 'Position', [100 100 1250 520]);
    ax = subplot(1,2,1);
    plot(ax, radialEccDegs, conesNumInPatch, 'rs-', 'MarkerSize', 10, 'MarkerFaceColor', [1 0.5 0.5], 'LineWidth', 1.5); 
    set(ax, 'XLim', [-1 31], 'XTick', 0:2:30);
    if (strcmp(analyzedRetinaMeridian,'temporal meridian'))
        set(ax, 'XDir', 'reverse')
    end
    xtickangle(ax, 0);
    grid(ax, 'on');
    legend(ax,{analyzedEye}, 'Location', 'northeast')
    ylabel(ax,'# of cones/deg2')
    xlabel(ax, sprintf('eccentricity, %s (degs)', analyzedRetinaMeridian));
    set(ax, 'FontSize', 16)

    ax = subplot(1,2,2);
    hold(ax, 'on');
    if (~isempty(visualToAnatomicalRcRatioPrcTilesOverSubjects))
        if (strcmp(analyzedRetinaMeridian, 'nasal meridian'))
            disp('nasal')
            n1 = 14;
            patch1 = patchCoordsFromXYmeanYstd(radialEccDegs(1:n1), visualToAnatomicalRcRatioPrcTilesOverSubjects(1:2,1:n1));
            patch1A = patchCoordsFromXYmeanYstd(radialEccDegs(1:n1), visualToAnatomicalRcRatioPrcTilesOverSubjects(3:4,1:n1));

            n2 = 19;
            patch2 = patchCoordsFromXYmeanYstd(radialEccDegs(n2:end), visualToAnatomicalRcRatioPrcTilesOverSubjects(1:2,n2:end));
            patch2A = patchCoordsFromXYmeanYstd(radialEccDegs(n2:end), visualToAnatomicalRcRatioPrcTilesOverSubjects(3:4,n2:end));
        else
            disp('temporal')
            patch1 = patchCoordsFromXYmeanYstd(radialEccDegs, visualToAnatomicalRcRatioPrcTilesOverSubjects(1:2,:));
            patch1A = patchCoordsFromXYmeanYstd(radialEccDegs, visualToAnatomicalRcRatioPrcTilesOverSubjects(3:4,:));
            patch2 = [];
            patch2A = [];
        end

        patch(ax, patch1.x, patch1.y, 0*patch1.y, [1 0 0], 'FaceColor', [1 0.75 0.75], 'EdgeColor', 'none', 'FaceAlpha', 0.5);
        patch(ax, patch1A.x, patch1A.y, 0*patch1A.y, [1 0 0], 'FaceColor', [1 0.5 0.5], 'EdgeColor', 'none', 'FaceAlpha', 0.5);
        if (~isempty(patch2))
            patch(ax, patch2.x, patch2.y, 0*patch2.y, [1 0 0], 'FaceColor', [1 0.75 0.75], 'EdgeColor', 'none','FaceAlpha', 0.5);
            patch(ax, patch2A.x, patch2A.y, 0*patch2A.y, [1 0 0], 'FaceColor', [1 0.5 0.5], 'EdgeColor', 'none','FaceAlpha', 0.5);
        end
        
    end

    
    plot(ax, radialEccDegs, visualToAnatomicalRcRatio, 'ro-', 'MarkerSize', 10, 'LineWidth', 1.5, 'MarkerFaceColor', [1 0.5 0.5]); hold on
    if (~isempty(visualToAnatomicalRcRatioPrcTilesOverSubjects))
        plot(ax, radialEccDegs, visualToAnatomicalRcRatioPrcTilesOverSubjects, 'k-');
    end

    set(ax, 'XLim', [-1 31], 'XTick', 0:2:30, 'YLim', [0 20], 'YTick', 0:1:30);
    xtickangle(ax, 0);
    if (strcmp(analyzedRetinaMeridian, 'temporal meridian'))
        set(ax, 'XDir', 'reverse')
    end

    grid(ax, 'on');

    ylabel(ax,'visual cone Rc/anatomical cone Rc')
    xlabel(ax, sprintf('eccentricity, %s (degs)', analyzedRetinaMeridian));
    set(ax, 'FontSize', 16)
    drawnow;

    NicePlot.exportFigToPDF(pdfFilename, hFig, 300);
    close(hFig);
end

function thePatch = patchCoordsFromXYmeanYstd(x, prctiles)
    minPrcTile = prctiles(1,:);
    maxPrcTile = prctiles(2,:);
    thePatch.x(1) = x(1);
    thePatch.y(1) = minPrcTile(1);

    thePatch.x(1+(1:numel(x))) = x;
    thePatch.y(1+(1:numel(x))) = maxPrcTile;

    thePatch.x(numel(thePatch.x)+1) = thePatch.x(numel(thePatch.x));
    thePatch.y(numel(thePatch.y)+1) = minPrcTile(end);

    thePatch.x(numel(thePatch.x)+(1:numel(x))) = fliplr(x);
    thePatch.y(numel(thePatch.y)+(1:numel(x))) = fliplr(minPrcTile);
end
