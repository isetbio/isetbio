function convolveGaussianWithPSF()
    
    retinalPoolingRadii = logspace(log10(0.01), log10(0.6), 8);
    plotlabOBJ = setupPlotLab();
    
    visualizeAnalysis = ~true;
    
    imposedRefractionErrorDiopters = 0.0; 
    eccTested = -[0 1 2.5 5.0 7.5 10 15 20 25];
    for eccIndex = 1:numel(eccTested)
        ecc = eccTested(eccIndex);
        doIt(ecc, imposedRefractionErrorDiopters, retinalPoolingRadii,  visualizeAnalysis, plotlabOBJ);
    end
    
end

function doIt(cellEcc, imposedRefractionErrorDiopters, retinalPoolingRadii,  visualizeAnalysis, plotlabOBJ )

    dataFileName = sprintf('cellEcc_%2.1f_cellRefractionError_%2.2fD_VisualGain.mat', cellEcc(1), imposedRefractionErrorDiopters);
    
    
    subjectsNum = 10;
    cMap = brewermap(512, 'greys');
    
    retinalRadius = zeros(numel(retinalPoolingRadii), subjectsNum);
    visualRadius = retinalRadius;
    visualGain = visualRadius;
    
    for kSubject = 1:(10*3)
        
        subjectID = mod(kSubject-1,10)+1;
        eccQuadrant = floor((kSubject-1)/10)+1;
        switch (eccQuadrant)
            case 1
                eccXrange = cellEcc(1)*[1 1];
                eccYrange = 0*[1 1];
            case 2
                eccYrange = cellEcc(1)*[1 1];
                eccXrange = 0*[1 1];
            case 3
                eccYrange = -cellEcc(1)*[1 1];
                eccXrange = 0*[1 1];
        end
        deltaEcc = 1;
    
         % Compute the subject PSF at the desired eccentricity
        [hEcc, vEcc, thePSFs, thePSFsupportDegs] = CronerKaplanRGCModel.psfAtEccentricity(subjectID, ...
                imposedRefractionErrorDiopters, eccXrange, eccYrange, deltaEcc);
            
        % Make a large PSF (via zero padding to be used to convolve with
        % the largest stimuli)
        thePSFOriginal = squeeze(thePSFs(1, 1, 1,:,:));
        thePSFsupportDegsOriginal = thePSFsupportDegs;
        
        psfSize = size(thePSFOriginal,1);
        largePSFsize = 1201;
        theLargePSF = zeros(largePSFsize, largePSFsize);
        margin = 0.5*(largePSFsize-psfSize);
        theLargePSF(margin+(1:psfSize), margin+(1:psfSize)) = thePSFOriginal;
        
        % Make corresponding large PSF support
        dS = thePSFsupportDegsOriginal(2)-thePSFsupportDegsOriginal(1);
        theLargePSFsupportDegs = -(0.5*(largePSFsize-1)*dS):dS:(0.5*(largePSFsize-1)*dS);
        
        
        % Convolve different retinal pooling regions and compute the
        % visually-mapped pooling region
        for retinalRadiusIndex = 1:numel(retinalPoolingRadii)
            
            rfPoolingRadiusDegsInRetinalSpace = retinalPoolingRadii(retinalRadiusIndex);
            if (rfPoolingRadiusDegsInRetinalSpace < 0.2)
                thePSF = thePSFOriginal;
                thePSFsupportDegs = thePSFsupportDegsOriginal;
            else
                thePSF = theLargePSF;
                thePSFsupportDegs = theLargePSFsupportDegs;
            end

            [X,Y] = meshgrid(thePSFsupportDegs,thePSFsupportDegs);
            rfPoolingInRetinalSpace = exp(-(X/rfPoolingRadiusDegsInRetinalSpace).^2).*exp(-(Y/rfPoolingRadiusDegsInRetinalSpace).^2);
            rfPoolingInRetinalSpaceNorm = rfPoolingInRetinalSpace / max(rfPoolingInRetinalSpace(:));
            [ellipseInRetinalSpace, semiAxesRetinalSpace, noFitRetinalSpace] = ...
                fitEllipse(thePSFsupportDegs, rfPoolingInRetinalSpaceNorm);

            rfPoolingInVisualSpace = conv2(rfPoolingInRetinalSpace, thePSF, 'same');
            rfPoolingInVisualSpaceNorm = rfPoolingInVisualSpace / max(rfPoolingInVisualSpace(:));
            [ellipseInVisualSpace, semiAxesVisualSpace, noFitVisualSpace] = ...
                fitEllipse(thePSFsupportDegs, rfPoolingInVisualSpaceNorm);

            % Saved data
            retinalRadius(retinalRadiusIndex, kSubject) = mean(semiAxesRetinalSpace);
            visualRadius(retinalRadiusIndex, kSubject) = mean(semiAxesVisualSpace);
            visualGain(retinalRadiusIndex, kSubject) = max(rfPoolingInVisualSpace(:))/max(rfPoolingInRetinalSpace(:));
            
            if ( visualizeAnalysis)
                % Extract profile at peak
                row = round(size(rfPoolingInRetinalSpace,1)/2);
                rfPoolingRetinalSpaceProfile = squeeze(rfPoolingInRetinalSpace(row,:));
                rfPoolingVisualSpaceProfile = squeeze(rfPoolingInVisualSpace(row,:));
                
                visualizeFit(thePSFsupportDegs, thePSF, subjectID, eccXrange, eccYrange, ...
                    rfPoolingInRetinalSpaceNorm, rfPoolingInVisualSpaceNorm, ...
                    rfPoolingRetinalSpaceProfile, rfPoolingVisualSpaceProfile, ...
                    retinalRadius(retinalRadiusIndex, kSubject), visualRadius(retinalRadiusIndex, kSubject), ...
                    visualGain(retinalRadiusIndex, kSubject), ellipseInRetinalSpace, ellipseInVisualSpace, cMap, plotlabOBJ );
            end
        end
    end
    save(dataFileName, 'retinalPoolingRadii', 'retinalRadius', 'visualRadius', 'visualGain');
end

function visualizeFit(thePSFsupportDegs, thePSF, subjectID, eccXrange, eccYrange, ...
                rfPoolingInRetinalSpace, rfPoolingInVisualSpace, ...
                rfPoolingRetinalSpaceProfile, rfPoolingVisualSpaceProfile, ...
                retinalRadius, visualRadius, ...
                visualGain, ellipseInRetinalSpace, ellipseInVisualSpace, cMap, plotlabOBJ )
            
        % Visualize them
        hFig = figure(1); clf;

        subplotPosVectors = NicePlot.getSubPlotPosVectors(...
               'rowsNum', 2, ...
               'colsNum', 4, ...
               'heightMargin',  0.02, ...
               'widthMargin',    0.02, ...
               'leftMargin',     0.03, ...
               'rightMargin',    0.00, ...
               'bottomMargin',   0.04, ...
               'topMargin',      0.02);

        xLims = 0.2*[-1 1];
        % Plot the PSF
        
        ax = subplot('Position', subplotPosVectors(1,1).v);
        theTitle = sprintf('subject #%d PSF at %2.1f, %2.1f degs', subjectID, eccXrange(1), eccYrange(1));
        renderPSF(ax, thePSFsupportDegs, thePSF, xLims, cMap, theTitle)

        ax = subplot('Position', subplotPosVectors(1,2).v);
        renderKernel(ax, thePSFsupportDegs, rfPoolingInRetinalSpace, xLims, cMap, 'pooling in retinal space');
        
        ax = subplot('Position', subplotPosVectors(2,2).v);
        renderKernelProfile(ax, thePSFsupportDegs, rfPoolingRetinalSpaceProfile, [], retinalRadius, [], xLims);
        
        ax = subplot('Position', subplotPosVectors(1,3).v);
        renderKernel(ax, thePSFsupportDegs, rfPoolingInVisualSpace, xLims, cMap, 'pooling in visual space');
            
        ax = subplot('Position', subplotPosVectors(2,3).v);
        renderKernelProfile(ax, thePSFsupportDegs, [], rfPoolingVisualSpaceProfile, [], visualRadius, xLims);

        ax = subplot('Position', subplotPosVectors(1,4).v);
        renderFittedEllipses(ax, ellipseInRetinalSpace, ellipseInVisualSpace, xLims, 'fitted ellipses');
        
        ax = subplot('Position', subplotPosVectors(2,4).v);
        renderKernelProfile(ax, thePSFsupportDegs, rfPoolingRetinalSpaceProfile, rfPoolingVisualSpaceProfile, [], [], xLims);

        
        ax = subplot('Position', subplotPosVectors(2,1).v);
        bar(ax, 0, visualGain, 1, 'r');
        set(ax, 'XLim', [-0.1 0.1], 'YLim', [0 1], 'XTick', [0], 'YTick', 0:0.2:1);
        axis(ax, 'square');
        
        %plotlabOBJ.exportFig(hFig, 'pdf', sprintf('subj_%d_ecc_%2.0f_%2.0f_radius_%2.3fDegs', subjectID, eccXrange(1), eccYrange(1), rfPoolingInRetinalSpace), pwd());

end

function renderPSF(ax, thePSFsupportDegs, thePSF, xLims, cMap, theTitle)
        imagesc(ax, thePSFsupportDegs, thePSFsupportDegs, thePSF); hold(ax, 'on');
        plot(ax, [0 0], [-1 1], 'k-'); plot(ax, [-1 1], [0 0], 'k-');
        axis(ax, 'square');  axis(ax, 'xy');
        set(ax, 'XLim', xLims, 'YLim', xLims);
        set(ax, 'XTickLabel', {}, 'YTickLabel', {});
        colormap(ax,cMap);
        title(ax, theTitle);
end

function renderKernel(ax, thePSFsupportDegs, rfPoolingInRetinalSpace, xLims, cMap, theTitle)
        contourf(ax, thePSFsupportDegs, thePSFsupportDegs, rfPoolingInRetinalSpace, 6);  hold(ax, 'on');
        plot(ax, [0 0], [-1 1], 'k-'); plot(ax, [-1 1], [0 0], 'k-');
        axis(ax, 'square');  axis(ax, 'xy');
        set(ax, 'XLim', xLims, 'YLim', xLims, 'CLim', [0 1.1]);
        colormap(ax,cMap);
        set(ax, 'XTickLabel', {}, 'YTickLabel', {}, 'XTick', -1:0.05:1, 'YTick', -1:0.05:1);
        axis(ax, 'square');
        title(ax, theTitle);
end


function renderKernelProfile(ax, thePSFsupportDegs, rfPoolingRetinalSpaceProfile, rfPoolingVisualSpaceProfile, retinalRadius, visualRadius, xLims)
        hold(ax, 'on');
        if (~isempty(rfPoolingVisualSpaceProfile))
            area(ax, thePSFsupportDegs, rfPoolingVisualSpaceProfile, 'LineWidth', 1.0, 'FaceColor', [1 0.5 0.5]);
        end
        if (~isempty(rfPoolingRetinalSpaceProfile))
            area(ax, thePSFsupportDegs, rfPoolingRetinalSpaceProfile,  'LineWidth', 1.0, 'FaceColor', [0.5 0.5 0.5]);
        end
        if (~isempty(retinalRadius))
            plot(retinalRadius*[-1 1], max(rfPoolingRetinalSpaceProfile)*exp(-1)*[1 1], 'k-');
        end
        if (~isempty(visualRadius))
            [~,idx] = max(rfPoolingVisualSpaceProfile);
            xo = thePSFsupportDegs(idx);
            plot(visualRadius*[-1 1]+xo, max(rfPoolingVisualSpaceProfile)*exp(-1)*[1 1], 'k-');
        end
        
        plot(ax, [0 0], [0 1], 'k-');
        axis(ax, 'square'); axis(ax, 'xy');
        set(ax, 'XLim', xLims, 'YLim', [0 1]);
        set(ax, 'XTick', -1:0.02:1, 'XTick', -1:0.05:1, 'YTick', 0:0.2:1, 'YTickLabel', {});
        axis(ax, 'square');
end

function renderFittedEllipses(ax, ellipseInRetinalSpace, ellipseInVisualSpace, xLims, theTitle)
        plot(ax, ellipseInRetinalSpace.x, ellipseInRetinalSpace.y, 'k-', 'LineWidth', 1.5);
        hold(ax, 'on');
        plot(ax, ellipseInVisualSpace.x, ellipseInVisualSpace.y, 'r-', 'LineWidth', 1.5);
        plot(ax, [0 0], [-1 1], 'k-'); plot(ax, [-1 1], [0 0], 'k-');
        axis(ax, 'square');  axis(ax, 'xy');
        set(ax, 'XLim', xLims, 'YLim', xLims,  'XTick', -1:0.05:1, 'YTick', -1:0.05:1);
        set(ax, 'XTickLabel', {}, 'YTickLabel', {});
        title(ax, theTitle);
end

function [ellipse, semiAxes, noFit] = fitEllipse(support, kernel)

    C = contourc(support, support, kernel, exp(-1)*[1 0.1]);
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

function plotlabOBJ = setupPlotLab()
    plotlabOBJ = plotlab();
    plotlabOBJ.applyRecipe(...
            'colorOrder', [0 0 0], ...
            'lineColor', [0.5 0.5 0.5], ...
            'lineWidth', 1.0, ...
            'scatterMarkerEdgeColor', [0.3 0.3 0.9], ...
            'lineMarkerSize', 10, ...
            'axesBox', 'off', ...
            'axesTickDir', 'in', ...
            'renderer', 'opengl', ...
            'axesTickLength', [0.01 0.01], ...
            'legendLocation', 'SouthWest', ...
            'figureWidthInches', 24, ...
            'figureHeightInches', 12);
end
