function convolveGaussianWithPSF()
    
    retinalPoolingRadii = logspace(log10(0.01), log10(0.6), 8);

    visualizeAnalysis = ~true;
    imposedRefractionErrorDiopters = 0.0;  % 0 == the defocus as measure
    for ecc = 0:40
        doIt([ecc 0], imposedRefractionErrorDiopters, retinalPoolingRadii,  visualizeAnalysis);
    end
    
    imposedRefractionErrorDiopters = 0.25; 
    for ecc = 0:40
        doIt([ecc 0], imposedRefractionErrorDiopters, retinalPoolingRadii,  visualizeAnalysis);
    end
    
end

function doIt(cellEcc, imposedRefractionErrorDiopters, retinalPoolingRadii,  visualizeAnalysis)

    dataFileName = sprintf('cellEcc_%2.1f_cellRefractionError_%2.2fD.mat', cellEcc(1), imposedRefractionErrorDiopters);
    eccXrange = cellEcc(1)*[1 1];
    eccYrange = cellEcc(2)*[1 1];
    deltaEcc = 1;
    
    subjectsNum = 10;
    cMap = brewermap(512, 'greys');
    
    retinalRadius = zeros(numel(retinalPoolingRadii), subjectsNum);
    visualRadius = retinalRadius;
    
    for subjectID = 1:10
         % Compute the subject PSF at the desired eccentricity
        [hEcc, vEcc, thePSFs, thePSFsupportDegs] = CronerKaplanRGCModel.psfAtEccentricity(subjectID, ...
                imposedRefractionErrorDiopters, eccXrange, eccYrange, deltaEcc);
            
        % Convolve different retinal pooling regions and compute the
        % visually-mapped pooling region
        for retinalRadiusIndex = 1:numel(retinalPoolingRadii)
            
            rfPoolingRadiusDegsInRetinalSpace = retinalPoolingRadii(retinalRadiusIndex);

            [X,Y] = meshgrid(thePSFsupportDegs,thePSFsupportDegs);
            rfPoolingInRetinalSpace = exp(-(X/rfPoolingRadiusDegsInRetinalSpace).^2).*exp(-(Y/rfPoolingRadiusDegsInRetinalSpace).^2);
            rfPoolingInRetinalSpace = rfPoolingInRetinalSpace / max(rfPoolingInRetinalSpace(:));
            [ellipseInRetinalSpace, semiAxesRetinalSpace, noFitRetinalSpace] = fitEllipse(thePSFsupportDegs, rfPoolingInRetinalSpace);

            rfPoolingInVisualSpace = conv2(rfPoolingInRetinalSpace, squeeze(thePSFs(1, 1, 1,:,:)), 'same');
            rfPoolingInVisualSpace = rfPoolingInVisualSpace / max(rfPoolingInVisualSpace(:));
            [ellipseInVisualSpace, semiAxesVisualSpace, noFitVisualSpace] = fitEllipse(thePSFsupportDegs, rfPoolingInVisualSpace);

            retinalRadius(retinalRadiusIndex, subjectID) = mean(semiAxesRetinalSpace);
            visualRadius(retinalRadiusIndex, subjectID) = mean(semiAxesVisualSpace);

            if ( visualizeAnalysis)
                visualizeFit(thePSFsupportDegs, thePSFs, subjectID, eccXrange, eccYrange, ...
                rfPoolingInRetinalSpace, rfPoolingInVisualSpace, ...
                retinalRadius(retinalRadiusIndex, subjectID), visualRadius(retinalRadiusIndex, subjectID), ...
                ellipseInRetinalSpace, ellipseInVisualSpace, cMap);
            end
        end
    end
    save(dataFileName, 'retinalPoolingRadii', 'retinalRadius', 'visualRadius');
    
    
    figure(2);
    subplot(1,3,1)
    plot(retinalPoolingRadii, retinalRadius, 'ko'); 
    axis 'equal'; axis 'square';
    xlabel('retinal radius');
    ylabel('retinal radius');
    
    subplot(1,3,2)
    plot(retinalPoolingRadii, visualRadius, 'bo');
    axis 'equal'; axis 'square';
    xlabel('retinal radius');
    ylabel('visual radius');
    
    subplot(1,3,3)
    plot(retinalRadius, visualRadius, 'ko');
    axis 'equal'; axis 'square';
    xlabel('retinal radius');
    ylabel('visual radius');
end

function visualizeFit(thePSFsupportDegs, thePSFs, subjectID, eccXrange, eccYrange, ...
                rfPoolingInRetinalSpace, rfPoolingInVisualSpace, retinalRadius, visualRadius, ...
                ellipseInRetinalSpace, ellipseInVisualSpace, cMap)
            
        % Visualize them
        hFig = figure(1); clf;
        set(hFig, 'Position', [10 10 1200 400]);

        subplotPosVectors = NicePlot.getSubPlotPosVectors(...
               'rowsNum', 1, ...
               'colsNum', 4, ...
               'heightMargin',  0.001, ...
               'widthMargin',    0.005, ...
               'leftMargin',     0.02, ...
               'rightMargin',    0.00, ...
               'bottomMargin',   0.02, ...
               'topMargin',      0.01);

        ax = subplot('Position', subplotPosVectors(1,1).v);
        imagesc(ax, thePSFsupportDegs, thePSFsupportDegs, squeeze(thePSFs(1, 1, 1,:,:))); hold(ax, 'on');
        plot(ax, [0 0], [-1 1], 'r-'); plot(ax, [-1 1], [0 0], 'r-');
        set(ax, 'XLim', 0.5*max(thePSFsupportDegs)*[-1 1], 'YLim', 0.5*max(thePSFsupportDegs)*[-1 1]);
        axis(ax, 'square');  axis(ax, 'xy');
        set(ax, 'XTickLabel', {}, 'YTickLabel', {});
        colormap(ax,cMap);
        title(ax, sprintf('Subject #%d PSF at %2.1f, %2.1f degs', subjectID, eccXrange(1), eccYrange(1)));

        ax = subplot('Position', subplotPosVectors(1,2).v);
        contourf(ax, thePSFsupportDegs, thePSFsupportDegs, rfPoolingInRetinalSpace, 10);  hold(ax, 'on');
        plot(ax, [0 0], [-1 1], 'r-'); plot(ax, [-1 1], [0 0], 'r-');
        set(ax, 'XLim', 0.5*max(thePSFsupportDegs)*[-1 1], 'YLim', 0.5*max(thePSFsupportDegs)*[-1 1]);
        axis(ax, 'square');  axis(ax, 'xy');
        set(ax, 'XTickLabel', {}, 'YTickLabel', {});
        axis(ax, 'square');
        colormap(ax,cMap);
        title(ax, 'Pooling in retinal space');

        ax = subplot('Position', subplotPosVectors(1,3).v);
        contourf(ax, thePSFsupportDegs, thePSFsupportDegs, rfPoolingInVisualSpace, 10);
        hold(ax, 'on');
        plot(ax, [0 0], [-1 1], 'r-'); plot(ax, [-1 1], [0 0], 'r-');
        set(ax, 'XLim', 0.5*max(thePSFsupportDegs)*[-1 1], 'YLim', 0.5*max(thePSFsupportDegs)*[-1 1]);
        axis(ax, 'square'); axis(ax, 'xy');
        set(ax, 'XTickLabel', {}, 'YTickLabel', {});
        axis(ax, 'square');
        colormap(ax,cMap);
        title(ax, 'pooling in visual space');

        ax = subplot('Position', subplotPosVectors(1,4).v);
        plot(ax, ellipseInRetinalSpace.x, ellipseInRetinalSpace.y, 'k-', 'LineWidth', 1.5);
        hold(ax, 'on');
        plot(ax, ellipseInVisualSpace.x, ellipseInVisualSpace.y, 'b-', 'LineWidth', 1.5);
        plot(ax, [0 0], [-1 1], 'r-'); plot(ax, [-1 1], [0 0], 'r-');
        set(ax, 'XLim', 0.5*max(thePSFsupportDegs)*[-1 1], 'YLim', 0.5*max(thePSFsupportDegs)*[-1 1]);
        axis(ax, 'square');  axis(ax, 'xy');
        set(ax, 'XTickLabel', {}, 'YTickLabel', {});
        axis(ax, 'square');
        legend({sprintf('retinal space radius: %2.4f degs', retinalRadius), ...
                sprintf('visual space radius: %2.4f degs', visualRadius)});
        colormap(ax,cMap);
        title(ax, 'fitted ellipses');
        pause(0.1)
end

function [ellipse, semiAxes, noFit] = fitEllipse(support, kernel)
    C = contour(support, support, kernel, exp(-1)*[1 0.1]);
    kk = 1;
    contourIndex = 0;
    while (kk < size(C,2))
        contourIndex = contourIndex+1;
        theLevel = C(1,kk);
       
        verticesNumForThisLevel = C(2,kk);
         if (theLevel == exp(-1))
            cM.level = theLevel;
            cM.x = C(1, kk+(1:verticesNumForThisLevel));
            cM.y = C(2, kk+(1:verticesNumForThisLevel));
         end
         
        kk = kk + verticesNumForThisLevel+1;
    end
    
    [ellipse.x,  ellipse.y, semiAxes, rfCenter, noFit] = fitEllipseToContour(cM.x,  cM.y);
    
end
