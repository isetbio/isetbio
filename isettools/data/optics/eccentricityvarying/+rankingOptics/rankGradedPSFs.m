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

    for includedSubjectCount = 1:numel(subjectPSFData)
        d = subjectPSFData{includedSubjectCount};
        if (isempty(d))
            continue;
        end
        subjectIDs(numel(subjectIDs)+1) = d.subjectID;

        % PSF cutoff frequency 
        cutoffFrequencies(numel(cutoffFrequencies)+1) = sqrt(d.psfXCutoffSF * d.psfYCutoffSF);
    end

    % Rank according to PSF cutoff frequency
    [cutoffFrequencies,sortedSubjectIndices] = sort(cutoffFrequencies, 'descend');
    rankedSubjectIDs = subjectIDs(sortedSubjectIndices);

    hFig = figure(1); clf;
    set(hFig, 'Position', [10 10 2000 700], 'Color', [1 1 1]);
    ax = subplot('Position', [0.05 0.05 0.94 0.94]);
    plot(ax, 1:numel(rankedSubjectIDs), cutoffFrequencies, 'o', 'MarkerSize', 16, 'MarkerFaceColor', [1 0.5 0.5], 'MarkerEdgeColor', [1 0 0]);
    set(ax, 'XLim', [0 numel(rankedSubjectIDs)+1], 'XTick', 1:numel(rankedSubjectIDs), 'XTickLabel', sprintf('%2.0f\n', rankedSubjectIDs));
    set(ax, 'YLim', [0 80], 'YTick', 0:10:100);
    set(ax, 'FontSize', 20);
    grid(ax, 'on');
    xtickangle(ax, 0);
    xlabel(ax, sprintf('%s subject ID', analyzedOpticsParams.zernikeDataBase));
    ylabel(ax, 'PSF cutoff frequency (c/deg)');

    switch (analyzedOpticsParams.whichEye)
        case 'right eye'
            titleString = sprintf('ranking based on PSF from RE\nat (x,y) ecc (degs): %2.1f, %2.1f', xPos, yPos);
        case 'left eye'
            titleString = sprintf('ranking based on PSF from LE\nat (x,y) ecc (degs): %2.1f, %2.1f', xPos, yPos);
    end
    title(ax, titleString, 'FontWeight', 'Normal');

    if (numel(rankedSubjectIDs) <= 12)
        colsNum = 4; rowsNum = 3;
        heightMargin = 0.1;
        bottomMargin =  0.05;
        topMargin = 0.04;
    elseif (numel(rankedSubjectIDs) <= 20)
        colsNum = 5; rowsNum = 4;
        heightMargin = 0.08;
        bottomMargin =  0.04;
        topMargin = 0.04;
    elseif (numel(rankedSubjectIDs) <= 40)
        colsNum = 8; rowsNum = 5;
        heightMargin = 0.07;
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
        bottomMargin =  0.03;
        topMargin = 0.03;
    elseif (numel(rankedSubjectIDs) <= 120)
        colsNum = 12; rowsNum = 10;
        heightMargin = 0.04;
        bottomMargin =  0.02;
        topMargin = 0.02;
    elseif (numel(rankedSubjectIDs) <= 150)
        colsNum = 15; rowsNum = 10;
        heightMargin = 0.04;
        bottomMargin =  0.02;
        topMargin = 0.02;
    else
        colsNum = 20; rowsNum = 10;
        heightMargin = 0.03;
        bottomMargin =  0.02;
        topMargin = 0.01;
    end


    sv = NicePlot.getSubPlotPosVectors(...
           'colsNum', colsNum, ...
           'rowsNum', rowsNum, ...
           'heightMargin',  heightMargin, ...
           'widthMargin',    0.01, ...
           'leftMargin',     0.02, ...
           'rightMargin',    0.00, ...
           'bottomMargin',   bottomMargin, ...
           'topMargin',      topMargin); 

    hFig = figure(2); clf;
    set(hFig, 'Position', [10 10 1900 1150], 'Color', [1 1 1]);
    cLUT = 1-gray(1024);
    colormap(cLUT);

    zLevels = 0:0.2:1;
    xTicks = -10:2:10;
    yTicks = -10:2:10;
    xLim  = 6*[-1 1];
    yLim = 6*[-1 1];


    for includedSubjectCount = 1:numel(rankedSubjectIDs)
    	row = floor((includedSubjectCount-1)/colsNum) + 1;
    	col = mod((includedSubjectCount-1),colsNum)+1;
    	% Get subject PSF
        d = subjectPSFData{sortedSubjectIndices(includedSubjectCount)};

        ax = subplot('Position', sv(row,col).v);
        contourf(ax, d.psfSupportArcMin, d.psfSupportArcMin, d.psfImage/max(d.psfImage(:)), zLevels);
        hold(ax, 'on');
        plot(ax, d.psfSupportArcMin, d.psfSupportArcMin*0, 'r-', 'LineWidth', 1.0);
        plot(ax, d.psfSupportArcMin*0, d.psfSupportArcMin, 'r-', 'LineWidth', 1.0);
        hold(ax, 'off');
        colorbar(ax)
        axis(ax, 'image'); axis(ax, 'xy'); grid(ax, 'on');
        set(ax, 'CLim', [0 1], 'ZLim', [0 1], 'XLim', xLim, 'YLim', yLim, 'XTick', xTicks, 'YTick', yTicks);
        set(ax, 'FontSize', 16);
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
                titleString = sprintf('%s, #%d (RE) \n(cutoff: %2.1f c/deg)', ...
                analyzedOpticsParams.zernikeDataBase, rankedSubjectIDs(includedSubjectCount), cutoffFrequencies(includedSubjectCount));
            case 'left eye'
                titleString = sprintf('%s, #%d (LE) \n(cutoff: %2.1f c/deg)', ...
                analyzedOpticsParams.zernikeDataBase, rankedSubjectIDs(includedSubjectCount), cutoffFrequencies(includedSubjectCount));
        end

        title(ax, titleString, 'FontWeight', 'Normal');
        drawnow;
     end
end

