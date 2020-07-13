function performGaussianConvolutionWithPolansPSFanalysis(obj, deconvolutionOpticsParams, varargin)

    defaultEccTested = [0 0.25 0.5 1 1.5 2.0 2.5 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25];
    
    defaultRetinalPoolingRadii = logspace(log10(0.001), log10(0.6), 10);
    
    % Parse input
    p = inputParser;
    p.addParameter('eccTested', defaultEccTested);
    p.addParameter('retinalPoolingRadii', defaultRetinalPoolingRadii);
    p.parse(varargin{:});
            
    % Validate the deconvolutionOpticsParams
    obj.validateDeconvolutionOpticsParams(deconvolutionOpticsParams);
    
    eccTested = -fliplr(p.Results.eccTested);
    for eccIndex = 1:numel(eccTested)
        ecc = eccTested(eccIndex);
        doIt(obj.psfDeconvolutionDir, ecc, p.Results.retinalPoolingRadii, deconvolutionOpticsParams);
    end
end

function doIt(rootDir, cellEcc, retinalPoolingRadii, deconvolutionOpticsParams)

    % Whether to visualize the analysis
    visualizeAnalysis = false;
    subjectsNum = numel(deconvolutionOpticsParams.PolansWavefrontAberrationSubjectIDsToCompute);
    quadrantsNum = numel(deconvolutionOpticsParams.quadrantsToCompute);
    
    subjectIDs = deconvolutionOpticsParams.PolansWavefrontAberrationSubjectIDsToCompute;
    quadrants = deconvolutionOpticsParams.quadrantsToCompute;
    
    % Preallocate memory
    retinalRadius = zeros(numel(retinalPoolingRadii), quadrantsNum,  subjectsNum);
    visualRadius = retinalRadius;
    visualGain = visualRadius;
   
    deltaEcc = 1;
    for qIndex = 1:quadrantsNum
        eccQuadrant = quadrants{qIndex};
        switch (eccQuadrant)
            case 'horizontal'
                % horizontal meridian (nasal)
                eccXrange = cellEcc(1)*[1 1];
                eccYrange = 0*[1 1];
            case 'upper vertical'
                % vertical meridian (upper)
                eccYrange = cellEcc(1)*[1 1];
                eccXrange = 0*[1 1];
            case 'lower vertical'
                % vertical meridian (lower)
                eccYrange = -cellEcc(1)*[1 1];
                eccXrange = 0*[1 1];
            otherwise
                error('Unknown Polans quadrant: ''%s''.', eccQuadrant);
        end
        
        for sIndex = 1:subjectsNum
            subjectID = subjectIDs(sIndex);
            % Compute the subject PSF at the desired eccentricity
            
            if (cellEcc <= 15)
                wavefrontSpatialSamples = 701;
            else
                wavefrontSpatialSamples = 1001;
            end
            
            pupilDiameterMM = 3.0;
            wavelengthsListToCompute = [550];
            micronsPerDegree = []; % empty so as to compute for each eccentricity
            imposedRefractionErrorDiopters = 0;
            
            [hEcc, vEcc, thePSFs, thePSFsupportDegs] = CronerKaplanRGCModel.psfAtEccentricity(subjectID, ...
                imposedRefractionErrorDiopters, pupilDiameterMM, wavelengthsListToCompute, micronsPerDegree, ...
                wavefrontSpatialSamples, eccXrange, eccYrange, deltaEcc);
            fprintf('PSF computed');
            
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
            parfor retinalRadiusIndex = 1:numel(retinalPoolingRadii)
                fprintf('Computing spread for quadrant %s, radius %d/%d\n', eccQuadrant, retinalRadiusIndex, numel(retinalPoolingRadii));
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
                retinalRadius(retinalRadiusIndex,qIndex, sIndex) = mean(semiAxesRetinalSpace);
                visualRadius(retinalRadiusIndex,qIndex, sIndex) = mean(semiAxesVisualSpace);
                visualGain(retinalRadiusIndex,qIndex, sIndex) = max(rfPoolingInVisualSpace(:))/max(rfPoolingInRetinalSpace(:));
                
                if (visualizeAnalysis)
                    % Extract profile at peak
                    row = round(size(rfPoolingInRetinalSpace,1)/2);
                    rfPoolingRetinalSpaceProfile = squeeze(rfPoolingInRetinalSpace(row,:));
                    rfPoolingVisualSpaceProfile = squeeze(rfPoolingInVisualSpace(row,:));
                    maxX = 0.1;
                    visualizeFit(thePSFsupportDegs, thePSF, subjectID, eccXrange, eccYrange, ...
                        rfPoolingInRetinalSpaceNorm, rfPoolingInVisualSpaceNorm, ...
                        rfPoolingRetinalSpaceProfile, rfPoolingVisualSpaceProfile, ...
                        retinalRadius(retinalRadiusIndex,qIndex, sIndex), visualRadius(retinalRadiusIndex,qIndex, sIndex), ...
                        visualGain(retinalRadiusIndex,qIndex, sIndex), ellipseInRetinalSpace, ellipseInVisualSpace, cMap, maxX);
                end
            end
        end % sIndex
    end % qIndex
    
    % Save data
    dataFileName = fullfile(rootDir, sprintf('ecc_%2.1f_deconvolutions_refractionError_%2.2fD.mat', cellEcc(1), imposedRefractionErrorDiopters));
    fprintf('Saving data to %s\n', dataFileName);
    save(dataFileName, 'retinalPoolingRadii', 'retinalRadius', 'visualRadius', 'visualGain', 'subjectIDs', 'quadrants');   
end

function visualizeFit(thePSFsupportDegs, thePSF, subjectID, eccXrange, eccYrange, ...
                rfPoolingInRetinalSpace, rfPoolingInVisualSpace, ...
                rfPoolingRetinalSpaceProfile, rfPoolingVisualSpaceProfile, ...
                retinalRadius, visualRadius, ...
                visualGain, ellipseInRetinalSpace, ellipseInVisualSpace, cMap, maxX)
            
        % Visualize them
        hFig = figure(1); clf;

        subplotPosVectors = NicePlot.getSubPlotPosVectors(...
               'rowsNum', 2, ...
               'colsNum', 4, ...
               'heightMargin',  0.03, ...
               'widthMargin',    0.02, ...
               'leftMargin',     0.03, ...
               'rightMargin',    0.00, ...
               'bottomMargin',   0.04, ...
               'topMargin',      0.02);

        xLims = maxX*[-1 1];
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
        drawnow;
        
end

function renderPSF(ax, thePSFSupportDegs, thePSF, xLims, cMap, theTitle)
    peakV = max(thePSF(:));
    imagesc(ax, thePSFSupportDegs, thePSFSupportDegs, thePSF); hold(ax, 'on');
    zLevels = logspace(log10(0.05), log10(1), 6)*peakV;
    [X,Y] = meshgrid(thePSFSupportDegs, thePSFSupportDegs);
    contour(ax, X,Y, thePSF, zLevels, 'LineColor', [0.5 0.5 0.5]);
    midRow = round(size(thePSF,1)/2)+1;
    profile = thePSF(midRow,:);
    deltaX = xLims(2)-xLims(1);
    area(ax,thePSFSupportDegs, xLims(1) + deltaX*profile/peakV, xLims(1), 'FaceColor', [1 0.5 0.5], 'EdgeColor', [1 0 0]);
    %plot(ax, [0 0], [-1 1], 'k-'); plot(ax, [-1 1], [0 0], 'k-');
    axis(ax, 'square');  axis(ax, 'xy');
    set(ax, 'XLim', xLims, 'YLim', xLims, 'CLim', [0 peakV]);
    set(ax, 'XTickLabel', {}, 'YTickLabel', {});
    set(ax, 'XTick', -1:0.05:1, 'YTick', -1:0.05:1, 'XTickLabel', {}, 'YTickLabel', {});
    colormap(ax,cMap);
    title(ax, theTitle);
end

function renderKernel(ax, thePSFsupportDegs, rfPoolingInRetinalSpace, xLims, cMap, theTitle)
        contourf(ax, thePSFsupportDegs, thePSFsupportDegs, rfPoolingInRetinalSpace, 6);  hold(ax, 'on');
        plot(ax, [0 0], [-1 1], 'k-'); plot(ax, [-1 1], [0 0], 'k-');
        axis(ax, 'square');  axis(ax, 'xy');
        set(ax, 'XLim', xLims, 'YLim', xLims, 'CLim', [0 1.1]);
        colormap(ax,cMap);
        ticks = -1:0.05:1;
        set(ax, 'XTickLabel', {}, 'YTickLabel', {}, ...
            'XTick', ticks, 'YTick', ticks, ...
            'XTickLabel', ticks*60, 'YTick', ticks*60);
        axis(ax, 'square');
        title(ax, theTitla@Qe);
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
        ticks = -1:0.05:1;
        set(ax, 'XTick', ticks,  'XTickLabel', ticks*60, 'YTick', 0:0.2:1, 'YTickLabel', {});
        axis(ax, 'square');
end

function renderFittedEllipses(ax, ellipseInRetinalSpace, ellipseInVisualSpace, xLims, theTitle)
        plot(ax, ellipseInRetinalSpace.x, ellipseInRetinalSpace.y, 'k-', 'LineWidth', 1.5);
        hold(ax, 'on');
        plot(ax, ellipseInVisualSpace.x, ellipseInVisualSpace.y, 'r-', 'LineWidth', 1.5);
        plot(ax, [0 0], [-1 1], 'k-'); plot(ax, [-1 1], [0 0], 'k-');
        axis(ax, 'square');  axis(ax, 'xy');
        ticks = -1:0.05:1;
        set(ax, 'XLim', xLims, 'YLim', xLims,  'XTick', ticks, 'YTick', ticks, ...
            'XTickLabel', ticks*60, 'YTickLabel', ticks*60);
        set(ax, 'XTickLabel', {}, 'YTickLabel', {});
        title(ax, theTitle);
end

function [ellipse, semiAxes, noFit] = fitEllipse(support, kernel)

    [~, idx] = max(kernel(:));
    [X,Y] = meshgrid(support, support);
    X = X(:);
    Y = Y(:);
    peakPos = [X(idx) Y(idx)];

    
    zLevels = exp(-1)*[1 0.9];
    C = contourc(support, support, kernel, zLevels);
    kk = 1;
    contourIndex = 0;
    centerNotEncountered = true;
    while (kk < size(C,2))
        contourIndex = contourIndex+1;
        theLevel = C(1,kk);
       
        verticesNumForThisLevel = C(2,kk);
        if (theLevel == zLevels(1)) && (centerNotEncountered)
            cM.level = theLevel;
            cM.x = C(1, kk+(1:verticesNumForThisLevel));
            cM.y = C(2, kk+(1:verticesNumForThisLevel));
        
            if (inpolygon(peakPos(1), peakPos(2), cM.x, cM.y))
                centerNotEncountered = false;
            end
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
            'renderer', 'painters', ...
            'axesTickLength', [0.01 0.01], ...
            'legendLocation', 'SouthWest', ...
            'figureWidthInches', 24, ...
            'figureHeightInches', 12);
end