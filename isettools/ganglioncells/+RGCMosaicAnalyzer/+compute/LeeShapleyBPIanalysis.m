%
% RGCMosaicAnalyzer.compute.LeeShapleyBPIanalysis()
%

function LeeShapleyBPIanalysis(...
	theMRGCMosaicSTFResponsesFullFileName,...
    theAnalyzedSTFsFullFileName)

	if (~isempty(strfind(theAnalyzedSTFsFullFileName, 'LconeIsolating')))
		theLconeIsolatingAnalyzedSTFFilename = theAnalyzedSTFsFullFileName;
		theMconeIsolatingAnalyzedSTFFilename = strrep(theAnalyzedSTFsFullFileName, 'LconeIsolating', 'MconeIsolating');
		theAchromaticAnalyzedSTFFilename = strrep(theAnalyzedSTFsFullFileName, 'LconeIsolating', 'Achromatic');
	elseif (~isempty(strfind(theAnalyzedSTFsFullFileName, 'MconeIsolating')))
		theMconeIsolatingAnalyzedSTFFilename = theAnalyzedSTFsFullFileName;
		theLconeIsolatingAnalyzedSTFFilename = strrep(theAnalyzedSTFsFullFileName, 'MconeIsolating', 'LconeIsolating');
		theAchromaticAnalyzedSTFFilename = strrep(theAnalyzedSTFsFullFileName, 'MconeIsolating', 'Achromatic');
	elseif (~isempty(strfind(theAnalyzedSTFsFullFileName, 'Achromatic')))
		theAchromaticAnalyzedSTFFilename = theAnalyzedSTFsFullFileName;
		theLconeIsolatingAnalyzedSTFFilename = strrep(theAnalyzedSTFsFullFileName, 'Achromatic', 'LconeIsolating');
		theMconeIsolatingAnalyzedSTFFilename = strrep(theAnalyzedSTFsFullFileName, 'Achromatic', 'MconeIsolating');
	else
		error('Input STF file is not LconeIsolating, MconeIsolating or Achromatic')
	end

	% Load the achromatic STF analyses
	d = load(theAchromaticAnalyzedSTFFilename);
	theAchromaticData.coneExcitationsBasedBPIs = d.theConeExcitationsBasedBPIs;
	theAchromaticData.photocurrentsBasedBPIs = d.thePhotocurrentsBasedBPIs;

	% Load the L-cone isolating STF analyses
	d = load(theLconeIsolatingAnalyzedSTFFilename);
	theLconeIsolatingData.coneExcitationsBasedBPIs = d.theConeExcitationsBasedBPIs;
	theLconeIsolatingData.photocurrentsBasedBPIs = d.thePhotocurrentsBasedBPIs;

	% Load the M-cone isolating STF analyses
	d = load(theMconeIsolatingAnalyzedSTFFilename);
	theMconeIsolatingData.coneExcitationsBasedBPIs = d.theConeExcitationsBasedBPIs;
	theMconeIsolatingData.photocurrentsBasedBPIs = d.thePhotocurrentsBasedBPIs;

	LconeCenterIdx = find(d.theCenterConeDominances == cMosaic.LCONE_ID);
    MconeCenterIdx = find(d.theCenterConeDominances == cMosaic.MCONE_ID);
	orientationsNum = numel(d.stimParams.orientationDegs);

	randomOrientationsLcenterCells = randi(orientationsNum,[1 numel(LconeCenterIdx)]);
	randomOrientationsMcenterCells = randi(orientationsNum,[1 numel(MconeCenterIdx)]);

	hFig = figure(); clf;
    set(hFig, 'Position', [10 10 1800 1000], 'Color', [1 1 1])

	for iOri = 1:orientationsNum+1

		% The cone excitations based STFs
        ax = subplot(2,orientationsNum+1,iOri);
        hold(ax, 'on')

        if (iOri <= orientationsNum)
	        y = theLconeIsolatingData.coneExcitationsBasedBPIs(LconeCenterIdx, iOri);
	        x = theAchromaticData.coneExcitationsBasedBPIs(LconeCenterIdx, iOri);
	        p1 = scatter(ax, x,y, 144, 'Marker', 'o', ...
	                'MarkerFaceColor', [1 0.5 0.5], 'MarkerEdgeColor', [0.5 0. 0.], ...
	                'MarkerFaceAlpha', 0.3, 'MarkerEdgeAlpha', 0.25, 'LineWidth', 1.0);

	        y = theMconeIsolatingData.coneExcitationsBasedBPIs(MconeCenterIdx, iOri);
	        x = theAchromaticData.coneExcitationsBasedBPIs(MconeCenterIdx, iOri);
	        p2 = scatter(ax, x,y, 144, 'Marker', 'o', ...
	                'MarkerFaceColor', [0.4 1.0 0.4], 'MarkerEdgeColor', [0. 0.5 0.], ...
	                'MarkerFaceAlpha', 0.3, 'MarkerEdgeAlpha', 0.25, 'LineWidth', 1.0);
	    else
	        % Now with the random orientations
	        for iiOri = 1:orientationsNum
	        	ii = find(randomOrientationsLcenterCells == iiOri);
	        	y = theLconeIsolatingData.coneExcitationsBasedBPIs(LconeCenterIdx(ii), iiOri);
	        	x = theAchromaticData.coneExcitationsBasedBPIs(LconeCenterIdx(ii), iiOri);
	        	if (iiOri == 1)
		        	p1 = scatter(ax, x, y, 144, 'Marker', 'o', ...
		                'MarkerFaceColor', [1 0.5 0.5], 'MarkerEdgeColor', [0.5 0. 0.], ...
		                'MarkerFaceAlpha', 0.3, 'MarkerEdgeAlpha', 0.25, 'LineWidth', 1.0);
		        else
		        	scatter(ax, x, y, 144, 'Marker', 'o', ...
		                'MarkerFaceColor', [1 0.5 0.5], 'MarkerEdgeColor', [0.5 0. 0.], ...
		                'MarkerFaceAlpha', 0.3, 'MarkerEdgeAlpha', 0.25, 'LineWidth', 1.0);
		        end
	        end

	        for iiOri = 1:orientationsNum
	        	ii = find(randomOrientationsMcenterCells == iiOri);
	        	y = theMconeIsolatingData.coneExcitationsBasedBPIs(MconeCenterIdx(ii), iiOri);
	        	x = theAchromaticData.coneExcitationsBasedBPIs(MconeCenterIdx(ii), iiOri);
	        	if (iiOri == 1)
		        	p2 = scatter(ax, x,y, 144, 'Marker', 'o', ...
		                'MarkerFaceColor', [0.4 1.0 0.4], 'MarkerEdgeColor', [0. 0.5 0.], ...
		                'MarkerFaceAlpha', 0.3, 'MarkerEdgeAlpha', 0.25, 'LineWidth', 1.0);
		        else
		        	scatter(ax, x,y, 144, 'Marker', 'o', ...
		                'MarkerFaceColor', [0.4 1.0 0.4], 'MarkerEdgeColor', [0. 0.5 0.], ...
		                'MarkerFaceAlpha', 0.3, 'MarkerEdgeAlpha', 0.25, 'LineWidth', 1.0);
		        end
	        end
	    end

        if (iOri == 1)
        	ylabel(ax, 'BPI (cone isolating)');
        end

        legend(ax, [p1 p2], {'L-cone dominated', 'M-cone dominated'}, 'Location', 'SouthEast');
        xtickangle(ax, 0);
        axis(ax, 'square')
        set(ax, 'XLim', [0 1], 'YLim', [0 1], 'XTick', 0:0.1:1, 'YTick', 0:0.1:1)
        grid(ax, 'on')
        if (iOri <= orientationsNum)
        	title(ax, sprintf('cone excitation STFs (%d degs)', d.stimParams.orientationDegs(iOri)));
        else
        	title(ax, sprintf('cone excitation STFs (random orientation)'));
        end

        set(ax, 'FontSize', 20);


        % The photocurrent based STFs
        ax = subplot(2, orientationsNum+1, orientationsNum+1+iOri);
        hold(ax, 'on')

        if (iOri <= orientationsNum)
	        y = theLconeIsolatingData.photocurrentsBasedBPIs(LconeCenterIdx, iOri);
	        x = theAchromaticData.photocurrentsBasedBPIs(LconeCenterIdx, iOri);
	        p1 = scatter(ax, x,y, 144, 'Marker', 'o', ...
	                'MarkerFaceColor', [1 0.5 0.5], 'MarkerEdgeColor', [0.5 0. 0.], ...
	                'MarkerFaceAlpha', 0.3, 'MarkerEdgeAlpha', 0.25, 'LineWidth', 1.0);

	        y = theMconeIsolatingData.photocurrentsBasedBPIs(MconeCenterIdx, iOri);
	        x = theAchromaticData.photocurrentsBasedBPIs(MconeCenterIdx, iOri);
	        p2 = scatter(ax, x,y, 144, 'Marker', 'o', ...
	                'MarkerFaceColor', [0.4 1.0 0.4], 'MarkerEdgeColor', [0. 0.5 0.], ...
	                'MarkerFaceAlpha', 0.3, 'MarkerEdgeAlpha', 0.25, 'LineWidth', 1.0);
	    else
	    	% Now with the random orientations
	        for iiOri = 1:orientationsNum
	        	ii = find(randomOrientationsLcenterCells == iiOri);
	        	y = theLconeIsolatingData.photocurrentsBasedBPIs(LconeCenterIdx(ii), iiOri);
	        	x = theAchromaticData.photocurrentsBasedBPIs(LconeCenterIdx(ii), iiOri);
	        	if (iiOri == 1)
		        	p1 = scatter(ax, x, y, 144, 'Marker', 'o', ...
		                'MarkerFaceColor', [1 0.5 0.5], 'MarkerEdgeColor', [0.5 0. 0.], ...
		                'MarkerFaceAlpha', 0.3, 'MarkerEdgeAlpha', 0.25, 'LineWidth', 1.0);
	        	else
	        		scatter(ax, x, y, 144, 'Marker', 'o', ...
		                'MarkerFaceColor', [1 0.5 0.5], 'MarkerEdgeColor', [0.5 0. 0.], ...
		                'MarkerFaceAlpha', 0.3, 'MarkerEdgeAlpha', 0.25, 'LineWidth', 1.0);
	        	end
	        end

	        for iiOri = 1:orientationsNum
	        	ii = find(randomOrientationsMcenterCells == iiOri);
	        	y = theMconeIsolatingData.photocurrentsBasedBPIs(MconeCenterIdx(ii), iiOri);
	        	x = theAchromaticData.photocurrentsBasedBPIs(MconeCenterIdx(ii), iiOri);
	        	if (iiOri == 1)
		        	p2 = scatter(ax, x,y, 144, 'Marker', 'o', ...
		                'MarkerFaceColor', [0.4 1.0 0.4], 'MarkerEdgeColor', [0. 0.5 0.], ...
		                'MarkerFaceAlpha', 0.3, 'MarkerEdgeAlpha', 0.25, 'LineWidth', 1.0);
		        else
		        	scatter(ax, x,y, 144, 'Marker', 'o', ...
		                'MarkerFaceColor', [0.4 1.0 0.4], 'MarkerEdgeColor', [0. 0.5 0.], ...
		                'MarkerFaceAlpha', 0.3, 'MarkerEdgeAlpha', 0.25, 'LineWidth', 1.0);
		        end
	        end
	    end


        xlabel(ax, 'BPI (achromatic)');
        legend(ax, [p1 p2], {'L-cone dominated', 'M-cone dominated'}, 'Location', 'SouthEast');
        xtickangle(ax, 0);
        axis(ax, 'square')
        set(ax, 'XLim', [0 1], 'YLim', [0 1], 'XTick', 0:0.1:1, 'YTick', 0:0.1:1)
        grid(ax, 'on')
        if (iOri <= orientationsNum)
        	title(ax, sprintf('photocurrent STFs (%d degs)', d.stimParams.orientationDegs(iOri)));
        else
        	title(ax, sprintf('photocurrent STFs (random orientation)'));
        end
        set(ax, 'FontSize', 20);


    end % for iORI

end
