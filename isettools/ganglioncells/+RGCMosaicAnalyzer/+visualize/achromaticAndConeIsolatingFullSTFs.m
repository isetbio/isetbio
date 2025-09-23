function achromaticAndConeIsolatingFullSTFs(figNo, theRGCindex, rfCenterConeDominance, optimalOrientation, sfSupport, orientationSupport, ...
          achromaticFullSTF, lConeIsolatingFullSTF, mConeIsolatingFullSTF, ...
          pdfFileName)

    maxSTFs = [max(achromaticFullSTF(:)) max(lConeIsolatingFullSTF(:)) max(mConeIsolatingFullSTF(:))]
    maxSTF = max(maxSTFs);

    renderFullSTF(figNo+1, rfCenterConeDominance, optimalOrientation, sfSupport, orientationSupport, achromaticFullSTF, maxSTF, ...
         sprintf('RGC #%d, achromatic', theRGCindex), pdfFileName);

    renderFullSTF(figNo+2,rfCenterConeDominance, optimalOrientation, sfSupport, orientationSupport, lConeIsolatingFullSTF, maxSTF, ...
         sprintf('RGC #%d, L-cone isolating', theRGCindex), pdfFileName);

    renderFullSTF(figNo+3, rfCenterConeDominance, optimalOrientation, sfSupport, orientationSupport, mConeIsolatingFullSTF, maxSTF, ...
         sprintf('RGC #%d, M-cone isolating', theRGCindex), pdfFileName);
	
end

function renderFullSTF(figNo, rfCenterConeDominance, optimalOrientation, sfSupport, orientationSupport, fullSTF, maxSTF, ...
    plotTitle, pdfFileName)

	ff = PublicationReadyPlotLib.figureComponents('1x1 standard figure');
	hFig = figure(figNo); clf;
    theAxes = PublicationReadyPlotLib.generatePanelAxes(hFig,ff);
    ax = theAxes{1,1};

    orientationSupport(end+1) = orientationSupport(1) + 180;
    fullSTF(end+1,:) = fullSTF(1,:);

    x = [];
    y = [];
    z = [];
    minSFincluded = 0.3;
    sfIndicesVisualized = [1 find(sfSupport>=minSFincluded)];
    sfSupport = sfSupport(sfIndicesVisualized);
    sfSupport(1) = 0;

    fullSTF = fullSTF(:,sfIndicesVisualized);

    for iSF = 1:numel(sfSupport)
        if (iSF == 1)
            x(1) = 0;
            y(1) = 0;
            z(1) = mean(squeeze(fullSTF(:, iSF)));
            xCPD(1) = sfSupport(iSF);
            yCPD(1) = sfSupport(iSF);
        else
            for iOri = 1:numel(orientationSupport)
                % Measured 
                x(numel(x)+1) = cosd(orientationSupport(iOri)) * (iSF-1);
                y(numel(y)+1) = sind(orientationSupport(iOri)) * (iSF-1);
                z(numel(z)+1) = fullSTF(iOri, iSF);
                xCPD(numel(xCPD)+1) = cosd(orientationSupport(iOri)) * sfSupport(iSF);
                yCPD(numel(yCPD)+1) = sind(orientationSupport(iOri)) * sfSupport(iSF);

                % Add ori + 180
                x(numel(x)+1) = cosd(orientationSupport(iOri)+180) * (iSF-1);
                y(numel(y)+1) = sind(orientationSupport(iOri)+180) * (iSF-1);
                z(numel(z)+1) = fullSTF(iOri, iSF);
                xCPD(numel(xCPD)+1) = cosd(orientationSupport(iOri)+180) * sfSupport(iSF);
                yCPD(numel(yCPD)+1) = sind(orientationSupport(iOri)+180) * sfSupport(iSF);

                % Negative frequency
                x(numel(x)+1) = -cosd(orientationSupport(iOri)) * (iSF-1);
                y(numel(y)+1) = -sind(orientationSupport(iOri)) * (iSF-1);
                z(numel(z)+1) = fullSTF(iOri, iSF);
                xCPD(numel(xCPD)+1) = -cosd(orientationSupport(iOri)) * sfSupport(iSF);
                yCPD(numel(yCPD)+1) = -sind(orientationSupport(iOri)) * sfSupport(iSF);

                x(numel(x)+1) = -cosd(orientationSupport(iOri)+180) * (iSF-1);
                y(numel(y)+1) = -sind(orientationSupport(iOri)+180) * (iSF-1);
                z(numel(z)+1) = fullSTF(iOri, iSF);
                xCPD(numel(xCPD)+1) = -cosd(orientationSupport(iOri)+180) * sfSupport(iSF);
                yCPD(numel(yCPD)+1) = -sind(orientationSupport(iOri)+180) * sfSupport(iSF);
            end % iOri
        end
    end % iSF

    % Make it ratio with respect to its amplitude at the lowest SF
    z = z / z(1);

    interpolationMethod = 'linear';
    extrapolationMethod = 'none';
    F = scatteredInterpolant(x(:),y(:),z(:),  interpolationMethod, extrapolationMethod);

    maxX = max(abs(x(:)));
    xx = -maxX:0.2:maxX;

    [xq,yq] = meshgrid(xx,xx);
    fullSTF2D = F(xq(:),yq(:));
    fullSTF2D = reshape(fullSTF2D, [numel(xx) numel(xx)]);

    zLevels = 0:0.1:2.0;

    contourf(ax, xx, xx, fullSTF2D, zLevels);
    hold(ax, 'on')
    
    xOptimalOrientation = maxX * cosd(optimalOrientation) * [-1 1];
    yOptimalOrientation = maxX * sind(optimalOrientation) * [-1 1];
    plot(ax, xOptimalOrientation, yOptimalOrientation, 'k-', 'LineWidth', 3.0);
    plot(ax, xOptimalOrientation, yOptimalOrientation, 'y--', 'LineWidth', 1.5);



    XLims = [xx(1)-1 xx(end)+1];

    sfsLabeled = [-60 -10 -3 0 3 10 60];
    for idx = 1:numel(sfsLabeled)
        [~,iSF] = min(abs(abs(sfsLabeled(idx))-sfSupport));
        if (sfsLabeled(idx)>0)
            if (abs(sfSupport(iSF))>9)
                XTickLabels{idx} = sprintf('%2.0f',round(sfSupport(iSF)/5)*5);
            else
                XTickLabels{idx} = sprintf('%2.0f', round(sfSupport(iSF)));
            end
            XTicks(idx) = iSF;
        elseif (sfsLabeled(idx)<0)
            if (abs(sfSupport(iSF))>9)
                XTickLabels{idx} = sprintf('%2.0f',-round(sfSupport(iSF)/5)*5);
            else
                XTickLabels{idx} = sprintf('%2.0f', round(sfSupport(iSF)));
            end
            XTicks(idx) = -iSF;
        elseif (sfsLabeled(idx) == 0)
            XTickLabels{idx} = '0';
            XTicks(idx) = 0;
        end
    end

    cMap = cat(1, brewermap(512, '*greys'), brewermap(512, 'reds'));

    axis(ax, 'image'); axis(ax, 'xy');
    set(ax, 'XLim', XLims, ...
            'XTick', XTicks, ...
            'XTickLabel', XTickLabels, ...
            'YLim', XLims, 'YTick', XTicks, ...
            'YTickLabel', XTickLabels, ...
            'Color', cMap(1,:), ...
            'CLim', [0 2], 'ZLim', [0 2]);

    xlabel(ax,'spatial frequency, x (c/deg)');
    ylabel(ax,'spatial frequency, y (c/deg)');



    switch (rfCenterConeDominance)
        case cMosaic.LCONE_ID
            theTitle = sprintf('%s (L-cone center)', plotTitle);
        case cMosaic.MCONE_ID
            theTitle =  sprintf('%s (M-cone center)', plotTitle);
        case cMosaic.SCONE_ID
            theTitle = sprintf('%s (S-cone center)', plotTitle);
    end
    title(ax, theTitle);
    colormap(ax,cMap);
    colorbar('Orientation', 'Vertical', 'Location', 'EastOutside')
    
    % Finalize figure using the Publication-Ready format
    PublicationReadyPlotLib.applyFormat(ax,ff);

    theRawFiguresDir = RGCMosaicConstructor.filepathFor.rawFigurePDFsDir();
    pdfExportSubDir = 'validation';
    thePDFfileName = fullfile(theRawFiguresDir, pdfExportSubDir, sprintf('%s_%s.pdf', pdfFileName, plotTitle));
    NicePlot.exportFigToPDF(thePDFfileName,hFig,  300);
end