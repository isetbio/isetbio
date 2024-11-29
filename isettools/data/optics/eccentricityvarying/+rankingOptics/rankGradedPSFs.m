function rankedSubjectIDs = rankGradedPSFs(exportsDir, analyzedOpticsParams)
	xPos = analyzedOpticsParams.eccentricityPosDegs.x;
	yPos = analyzedOpticsParams.eccentricityPosDegs.y;
	switch (analyzedOpticsParams.whichEye)
        case 'right eye'
           dataFile = sprintf('%sOpticsRE_AnalysisEccPosDegs_%2.1f_%2.1f.mat', ...
                analyzedOpticsParams.zernikeDataBase, xPos, yPos);
           case 'left eye'
           dataFile = sprintf('%sOpticsLE_AnalysisEccPosDegs_%2.1f_%2.1f.mat', ...
                analyzedOpticsParams.zernikeDataBase,  xPos, yPos);
    end
    load(fullfile(exportsDir, dataFile), 'subjectPSFData');

	subjectIDs = [];
    cutoffFrequencies = [];
    MTFs = [];

    for includedSubjectCount = 1:numel(subjectPSFData)
        d = subjectPSFData{includedSubjectCount};
        if (isempty(d))
            continue;
        end
        subjectIDs(numel(subjectIDs)+1) = d.subjectID;

        % PSF cutoff frequency 
        cutoffFrequencies(numel(cutoffFrequencies)+1) = sqrt(d.psfXCutoffSF * d.psfYCutoffSF);

        % The MTF spectra along the major and minor axes
        % Major axis
        MTFs(numel(cutoffFrequencies)+1,1,:) = d.spectralSupportCyclesPerDegreeX;
        MTFs(numel(cutoffFrequencies)+1,2,:) = d.PSFSpectrumSliceX;
        % Minor axis
        MTFs(numel(cutoffFrequencies)+1,3,:) = d.spectralSupportCyclesPerDegreeY;
        MTFs(numel(cutoffFrequencies)+1,4,:) = d.PSFSpectrumSliceY;
    end

    % Rank according to PSF cutoff frequency
    [cutoffFrequencies,sortedSubjectIndices] = sort(cutoffFrequencies, 'descend');
    rankedSubjectIDs = subjectIDs(sortedSubjectIndices);

    hFig = figure(1); clf;
    set(hFig, 'Position', [10 10 650 700], 'Color', [1 1 1]);
    ax = subplot('Position', [0.1 0.1 0.87 0.85]);
    plot(ax, 1:numel(rankedSubjectIDs), cutoffFrequencies, 'k-',  'LineWidth', 4.0);
    hold(ax, 'on');
    plot(ax, 1:numel(rankedSubjectIDs), cutoffFrequencies, 'o-', ...
        'MarkerFaceColor', [0.95 0.75 0.0], 'MarkerEdgeColor', [0.2 0.2 0.2], 'MarkerSize', 20, 'Color', [0.95 0.75 0.0], 'LineWidth', 2.0);
    set(ax, 'XLim', [0 numel(rankedSubjectIDs)+1], 'XTick', 1:numel(rankedSubjectIDs), 'XTickLabel', sprintf('%2.0f\n', rankedSubjectIDs));
    set(ax, 'YLim', [0 65], 'YTick', 0:10:100);
    set(ax, 'FontSize', 20);
    grid(ax, 'on');
    xtickangle(ax, 0);
    xlabel(ax, sprintf('subject ID (%s)', analyzedOpticsParams.zernikeDataBase));
    ylabel(ax, 'PSF corner (-15dB) frequency (c/deg)');

    switch (analyzedOpticsParams.whichEye)
        case 'right eye'
            titleString = sprintf('ranking based on PSF from RE at (x,y) ecc (degs): %2.1f, %2.1f', xPos, yPos);
        case 'left eye'
            titleString = sprintf('ranking based on PSF from LE at (x,y) ecc (degs): %2.1f, %2.1f', xPos, yPos);
    end
    title(ax, titleString, 'FontWeight', 'Normal');

    if (numel(rankedSubjectIDs) <= 12)
        colsNum = 4; rowsNum = 3;
        heightMargin = 0.1;
        widthMargin = 0.05;
        bottomMargin =  0.065;
        topMargin = 0.03;
    elseif (numel(rankedSubjectIDs) <= 20)
        colsNum = 5; rowsNum = 4;
        heightMargin = 0.08;
        widthMargin = 0.04;
        bottomMargin =  0.04;
        topMargin = 0.04;
    elseif (numel(rankedSubjectIDs) <= 40)
        colsNum = 8; rowsNum = 5;
        heightMargin = 0.07;
        widthMargin = 0.035;
        bottomMargin =  0.04;
        topMargin = 0.04;
    elseif (numel(rankedSubjectIDs) <= 60)
        colsNum = 10; rowsNum = 6;
        heightMargin = 0.06;
        bottomMargin =  0.03;
        topMargin = 0.04;
    elseif (numel(rankedSubjectIDs) <= 80)
        colsNum = 10; rowsNum = 8;
        heightMargin = 0.05;
        widthMargin = 0.03;
        bottomMargin =  0.03;
        topMargin = 0.03;
    elseif (numel(rankedSubjectIDs) <= 120)
        colsNum = 12; rowsNum = 10;
        heightMargin = 0.04;
        widthMargin = 0.025;
        bottomMargin =  0.02;
        topMargin = 0.02;
    elseif (numel(rankedSubjectIDs) <= 150)
        colsNum = 15; rowsNum = 10;
        heightMargin = 0.04;
        widthMargin = 0.025;
        bottomMargin =  0.02;
        topMargin = 0.02;
    else
        colsNum = 20; rowsNum = 10;
        heightMargin = 0.03;
        widthMargin = 0.02;
        bottomMargin =  0.02;
        topMargin = 0.01;
    end


    sv = NicePlot.getSubPlotPosVectors(...
           'colsNum', colsNum, ...
           'rowsNum', rowsNum, ...
           'heightMargin',  heightMargin, ...
           'widthMargin',    widthMargin, ...
           'leftMargin',     0.04, ...
           'rightMargin',    0.01, ...
           'bottomMargin',   bottomMargin, ...
           'topMargin',      topMargin); 

    
    hFig = figure(2); clf;
    set(hFig, 'Position', [10 10 1900 1150], 'Color', [1 1 1]);

    cornerFrequencyAttenuationDB = -15;
    cornerFrequencyAttenuation = 10^(cornerFrequencyAttenuationDB/20);

    for includedSubjectCount = 1:numel(rankedSubjectIDs)
        % Get subject PSF
        d = subjectPSFData{sortedSubjectIndices(includedSubjectCount)};

        row = floor((includedSubjectCount-1)/colsNum) + 1;
        col = mod((includedSubjectCount-1),colsNum)+1;
        ax = subplot('Position', sv(row,col).v);

        % Map the 0 c/deg to 0.1 c/deg
        minSF = 0.3;
        d.spectralSupportCyclesPerDegreeX(d.spectralSupportCyclesPerDegreeX==0) = minSF;
        d.spectralSupportCyclesPerDegreeY(d.spectralSupportCyclesPerDegreeY==0) = minSF;

        hold (ax, 'on');
        p1 = plot(ax, d.spectralSupportCyclesPerDegreeX, d.coneSpectrumSlice, 'k--', 'LineWidth', 2.0);
        plot(ax, d.spectralSupportCyclesPerDegreeX, d.coneSpectrumSlice, 'k-', 'LineWidth', 4.0);
        plot(ax, d.spectralSupportCyclesPerDegreeX, d.coneSpectrumSlice, 'w--', 'LineWidth', 2.0);
        plot(ax, d.spectralSupportCyclesPerDegreeX, d.PSFSpectrumSliceX, 'k-', 'LineWidth', 4.0);
        p2 = plot(ax, d.spectralSupportCyclesPerDegreeX, d.PSFSpectrumSliceX, 'r-', 'LineWidth', 2.0, 'Color', [1 0.3 0.3]);
        plot(ax, d.spectralSupportCyclesPerDegreeY, d.PSFSpectrumSliceY, 'k-', 'LineWidth', 4.0);
        p3 = plot(ax, d.spectralSupportCyclesPerDegreeY, d.PSFSpectrumSliceY, 'b-', 'LineWidth', 2.0, 'Color', [0.3 0.6 1.0]);
        plot(ax, d.spectralSupportCyclesPerDegreeX, sqrt(d.PSFSpectrumSliceX.*d.PSFSpectrumSliceY), 'k-', ...
            'LineWidth', 4.0);
        plot(ax, d.spectralSupportCyclesPerDegreeX, sqrt(d.PSFSpectrumSliceX.*d.PSFSpectrumSliceY), '-', ...
            'Color', [0.75 0.75 0.75], 'LineWidth', 2.0);
        plot(ax, d.spectralSupportCyclesPerDegreeX, cornerFrequencyAttenuation + 0*d.spectralSupportCyclesPerDegreeX, 'k:', 'LineWidth', 1.0);
        p4 = plot(ax, sqrt(d.psfXCutoffSF * d.psfYCutoffSF), cornerFrequencyAttenuation, 'ko', ...
            'MarkerFaceColor', [0.95 0.75 0.0], 'MarkerEdgeColor', [0.2 0.2 0.2], 'MarkerSize', 20, 'LineWidth', 1.5);
        if (includedSubjectCount == 1)
            legend(ax, [p1 p2 p3 p4], {'cones' 'cones + optics (minor)', 'cones + optics (major)', 'corner frequency (-15 dB)'}, ...
                'Location', 'SouthWest', 'FontSize', 14);
        end
        grid(ax, 'on');
        set(ax, 'XScale', 'log', 'XTick', [0.1 0.3 1 3 10 30 100 300], 'XLim', [minSF 300]);
        set(ax, 'YScale', 'linear', 'YTick', [0:0.2:1.0], 'YLim', [0 1]);
        set(ax, 'FontSize', 26);
        xtickangle(ax, 0);

        if (row == rowsNum)
            xlabel(ax, 'spatial frequency (c/deg)'); 
        end

        if (col == 1)
            ylabel(ax, 'power');
        else
            set(ax, 'YTickLabel', {});
        end

        switch (analyzedOpticsParams.whichEye)
            case 'right eye'
                titleString = sprintf('%s (RE) ID:#%d\ncornerF: %2.1f c/deg', ...
                analyzedOpticsParams.zernikeDataBase, rankedSubjectIDs(includedSubjectCount), cutoffFrequencies(includedSubjectCount));
            case 'left eye'
                titleString = sprintf('%s (LE) ID:#%d\ncornerF: %2.1f c/deg', ...
                analyzedOpticsParams.zernikeDataBase, rankedSubjectIDs(includedSubjectCount), cutoffFrequencies(includedSubjectCount));
        end

        title(ax, titleString, 'FontWeight', 'Normal', 'FontSize', 18);

    end % includedSubjectCount 


    hFig = figure(3); clf;
    set(hFig, 'Position', [10 10 1900 1150], 'Color', [1 1 1]);
    cLUT = 1-gray(1024);
    colormap(cLUT);

    zLevels = 0:0.2:1;
    xTicks = -10:2:10;
    yTicks = -10:2:10;
    xLim  = 6*[-1 1];
    yLim = 6*[-1 1];


    for includedSubjectCount = 1:numel(rankedSubjectIDs)
        % Get subject PSF
        d = subjectPSFData{sortedSubjectIndices(includedSubjectCount)};

    	row = floor((includedSubjectCount-1)/colsNum) + 1;
    	col = mod((includedSubjectCount-1),colsNum)+1;
        ax = subplot('Position', sv(row,col).v);
        contourf(ax, d.psfSupportArcMin, d.psfSupportArcMin, d.psfImage/max(d.psfImage(:)), zLevels);
        hold(ax, 'on');
        plot(ax, d.psfSupportArcMin, d.psfSupportArcMin*0, 'r-', 'LineWidth', 1.0);
        plot(ax, d.psfSupportArcMin*0, d.psfSupportArcMin, 'r-', 'LineWidth', 1.0);
        hold(ax, 'off');
        colorbar(ax)
        axis(ax, 'image'); axis(ax, 'xy'); grid(ax, 'on');
        set(ax, 'CLim', [0 1], 'ZLim', [0 1], 'XLim', xLim, 'YLim', yLim, 'XTick', xTicks, 'YTick', yTicks);
        set(ax, 'FontSize', 20);
        if (row == rowsNum)
        	xlabel(ax, 'arc min'); 
        else
        	set(ax, 'XTickLabel', {});
        end

        if (col == 1)
        	ylabel(ax, 'arc min');
        else
			set(ax, 'YTickLabel', {});
        end
        
        switch (analyzedOpticsParams.whichEye)
            case 'right eye'
                titleString = sprintf('%s (RE) ID:#%d', ...
                analyzedOpticsParams.zernikeDataBase, rankedSubjectIDs(includedSubjectCount));
            case 'left eye'
                titleString = sprintf('%s (LE) ID:#%d', ...
                analyzedOpticsParams.zernikeDataBase, rankedSubjectIDs(includedSubjectCount));
        end

        title(ax, titleString, 'FontWeight', 'Normal', 'FontSize', 18);
        drawnow;
     end
end

