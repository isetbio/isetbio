function batchGenerateRetinalParamsDataFilesForTargetVisualRF

    % Optics params
    ZernikeDataBase = 'Artal2012';
    examinedSubjectRankOrder = 1;
    pupilDiameterMM = 3.0;

    % Retinal location and eye
    analyzedRetinaMeridian = 'nasal meridian';
    
    % Number of cones in RF center
    conesNumPooledByTheRFcenter = 2;

    switch (conesNumPooledByTheRFcenter)
        case 1
            minEccForThisCenterConesNum = 0;
            maxEccForThisCenterConesNum = 11;
        case 2
            minEccForThisCenterConesNum = 0;
            maxEccForThisCenterConesNum = 11;
    end


    analyzedEye = 'right eye';
    subjectRankingEye = 'right eye';

    analysisFileName = sprintf('%s_subjRankOrder_%d_%s_%s_pupilDiamMM_%2.2f_conesNumPooledByRFcenter%d.mat', ...
        ZernikeDataBase, examinedSubjectRankOrder, ...
        strrep(analyzedEye, ' ', '_'), ...
        strrep(analyzedRetinaMeridian, ' ', '_'), ...
        pupilDiameterMM, conesNumPooledByTheRFcenter)
 
    regenerateData = true;
    if (regenerateData)

        % Sampling the eccentricity range. We sample very fine initially to
        % account for the fast reduction in cone density in the center. The
        % optics do not vary so fast.
        analyzedRadialEccDegs = [0];
        step = 0.1;
        while (analyzedRadialEccDegs(end) < 30)
            analyzedRadialEccDegs(numel(analyzedRadialEccDegs)+1) = ...
                analyzedRadialEccDegs(end) + step;
            step = step * 1.05;
            if (step>0.5)
                analyzedRadialEccDegs(end) = round(2*analyzedRadialEccDegs(end))/2;
                step = 0.5;
            end
        end

        analyzedRadialEccDegs = analyzedRadialEccDegs(...
            (analyzedRadialEccDegs >= minEccForThisCenterConesNum) & ...
            (analyzedRadialEccDegs <= maxEccForThisCenterConesNum)...
            );

        % From Croner & Kaplan '95 (Figure 4c and text)
        % "P surrounds were on average 6.7 times wider than the centers of 
        % the same cells, or about 45 times larger in area".
        surroundToCenterRcRatio = 6.7;
    
        % From Croner & Kaplan '95 (Figure 10b)
        % "These mean ratios for P and M cells are not significantly different 
        % (Student's t-test: P = 0.482). The overall mean ratio is 0.55. In 
        % general, then, the surround of a cell is less sensitive to contrast 
        % than is the center, and the ratio of sensitivities is such that the 
        % surround can reduce the center's response by about 55%.
        surroundToCenterIntegratedRatio = 0.55;

        % Make it large enough to encompass the largest PSFs in the periphery
        wavefrontSpatialSamples = 601;

        % Select between {'GaussianCenterDoubleExponentSurroundBased', 'GaussianCenterGaussianSurroundBased'}
        retinalConePoolingModel = 'GaussianCenterGaussianSurroundBased';

        % Struct with DoG model params for the target visual RF
        targetVisualRFDoGparams = struct(...
            'Kc', 1, ...
            'RcDegs', [], ...  % empty indicates that we should derive the center based on the conesNumPooledByTheRFcenter
            'conesNumPooledByTheRFcenter', conesNumPooledByTheRFcenter, ...
            'surroundToCenterRcRatio', surroundToCenterRcRatio, ...
            'surroundToCenterIntegratedRatio', surroundToCenterIntegratedRatio);
    
        % Struct with the various optics params
        opticsParams = struct(...
            'radialEccDegs', analyzedRadialEccDegs, ...
            'analyzedRetinaMeridian', analyzedRetinaMeridian, ...
            'ZernikeDataBase', ZernikeDataBase, ...
            'examinedSubjectRankOrder', examinedSubjectRankOrder , ...
            'analyzedEye', analyzedEye, ...
            'subjectRankingEye', subjectRankingEye, ...
            'pupilDiameterMM', pupilDiameterMM, ...
            'wavefrontSpatialSamples', wavefrontSpatialSamples ...
            );

        % Go !
        [retinalRFparamsDictionary, opticsParams, targetVisualRFDoGparams] = ...
            computeRetinalRFparamsAcrossEccentricities(opticsParams, targetVisualRFDoGparams, retinalConePoolingModel);
    
        % Save computed data
        save(analysisFileName, 'retinalRFparamsDictionary', 'opticsParams', 'targetVisualRFDoGparams', 'analyzedRadialEccDegs');
    else
        % Load computed data
        load(analysisFileName, 'retinalRFparamsDictionary', 'opticsParams', 'targetVisualRFDoGparams', 'analyzedRadialEccDegs');
    end

    % Evaluate generated RFs at a target eccentricity
    targetRadialEccDegs = 0.9;
    evaluteGeneratedRFs(retinalRFparamsDictionary, opticsParams, targetVisualRFDoGparams, targetRadialEccDegs);
end


function [retinalRFparamsDictionary, opticsParams, targetVisualRFDoGparams] = ...
    computeRetinalRFparamsAcrossEccentricities(opticsParams, targetVisualRFDoGparams, retinalConePoolingModel)

    % Instantiate the RetinaToVisualFieldTransformer
    xFormer = RetinaToVisualFieldTransformer('ZernikeDataBase', opticsParams.ZernikeDataBase);

    % Do the computation for all subjects and all eccentricities
    subjectRankOrder = opticsParams.examinedSubjectRankOrder;
    subjID = xFormer.subjectWithRankInEye(subjectRankOrder, opticsParams.subjectRankingEye);
    opticsParams.subjID = subjID;

    % Retrieve the horizontal and vertical eccs corresponding to the target
    % radial ecc, retinal meridian and eye
    [horizontalEccDegs, verticalEccDegs] = cMosaic.eccentricitiesForRetinaMeridianInEye(...
            opticsParams.radialEccDegs, opticsParams.analyzedRetinaMeridian, opticsParams.analyzedEye);

    retinalRFparamsDictionary = containers.Map();

    for iEcc = 1:numel(horizontalEccDegs)
        % Analyze effect of optics at this eccentricity
        eccDegs = [horizontalEccDegs(iEcc) verticalEccDegs(iEcc)];

        % Compute maxSpatialSupportDegs
        dStruct = xFormer.estimateConeCharacteristicRadiusInVisualSpace(...
                opticsParams.analyzedEye, eccDegs, subjID, opticsParams.pupilDiameterMM, '', ...
                'anatomicalConeCharacteristicRadiusDegs', [], ...
                'hFig', [], ...
                'videoOBJ', []);
        visualConeCharacteristicRadiusDegs = dStruct.visualConeCharacteristicRadiusDegs;
        maxSpatialSupportDegs = ...
            round((visualConeCharacteristicRadiusDegs * 2.0 * ...
                   targetVisualRFDoGparams.conesNumPooledByTheRFcenter * ...
                   targetVisualRFDoGparams.surroundToCenterRcRatio) * 100.0)/100;

        % Call retinalRFparamsForTargetVisualRF() to estimate etinalRFparamsStruct
        [retinalRFparamsStruct, weightsComputeFunctionHandle, ...
         targetVisualRF, spatialSupportDegs, theCircularPSFData] = xFormer.retinalRFparamsForTargetVisualRF(...
            targetVisualRFDoGparams, eccDegs, opticsParams.analyzedEye, subjID, opticsParams.pupilDiameterMM, ...
            'wavelengthSupportForVLambdaPSF', 550, ...
            'maxSpatialSupportDegs', maxSpatialSupportDegs, ...
            'wavefrontSpatialSamples', opticsParams.wavefrontSpatialSamples, ...
            'retinalConePoolingModel', retinalConePoolingModel, ... 
            'minimizationDomain', 'visual' ...  % Select between 'retinal' (involves deconvolution) and 'visual'
            );

        % Save to dictionary
        eccLabel = sprintf('Ecc_%2.3f_%2.3f', eccDegs(1), eccDegs(2));
        s.retinalConePoolingModel = retinalConePoolingModel;
        s.retinalRFparamsStruct = retinalRFparamsStruct;
        s.weightsComputeFunctionHandle = weightsComputeFunctionHandle;
        s.targetVisualRF = targetVisualRF;
        s.targetVisualRFspatialSupportDegs = spatialSupportDegs;
        s.maxSpatialSupportDegs = maxSpatialSupportDegs;
        s.theEmployedCircularPSFData = theCircularPSFData;
        retinalRFparamsDictionary(eccLabel) = s;
    end
end

function evaluteGeneratedRFs(retinalRFparamsDictionary, opticsParams, targetVisualRFDoGparams, targetRadialEccDegs)

    % Retrieve the horizontal and vertical eccs corresponding to the target
    % radial ecc, retinal meridian and eye
    [horizontalEccDegs, verticalEccDegs] = cMosaic.eccentricitiesForRetinaMeridianInEye(...
            opticsParams.radialEccDegs, opticsParams.analyzedRetinaMeridian, opticsParams.analyzedEye);

    [targetHorizontalEccDegs, targetVerticalEccDegs] = cMosaic.eccentricitiesForRetinaMeridianInEye(...
            targetRadialEccDegs, opticsParams.analyzedRetinaMeridian, opticsParams.analyzedEye);


    % Find the best source data set
    d = sqrt((horizontalEccDegs-targetHorizontalEccDegs(1)).^2+(verticalEccDegs-targetVerticalEccDegs(1)).^2);
    [~,iEcc] = min(d(:));
    sourceEccDegs = [horizontalEccDegs(iEcc) verticalEccDegs(iEcc)];

    fprintf('Will use the (%2.3f,%2.3f degs) dataset which is closest to the target radial eccentricity (%2.3f degs)\n', ...
        sourceEccDegs(1), sourceEccDegs(2), targetRadialEccDegs);

    % Retrieve computed data for the source ecc
    eccLabel = sprintf('Ecc_%2.3f_%2.3f', sourceEccDegs(1), sourceEccDegs(2));
    s = retinalRFparamsDictionary(eccLabel);
    retinalRFparamsStruct = s.retinalRFparamsStruct;
    weightsComputeFunctionHandle = s.weightsComputeFunctionHandle;
    targetVisualRF = s.targetVisualRF;
    rfSpatialSupportDegs = s.targetVisualRFspatialSupportDegs;
    maxSpatialSupportDegs = s.maxSpatialSupportDegs;
    theEmployedPSFData = s.theEmployedCircularPSFData;
    clear 's';

        
    % Generate a @cMosaic object located at the target eccentricity
    coneMosaicSize = max([0.5 2*maxSpatialSupportDegs]);
    targetConeMosaic = cMosaic(...
            'whichEye', opticsParams.analyzedEye, ...
            'sizeDegs', [1 1] * coneMosaicSize, ...
            'eccentricityDegs', [targetHorizontalEccDegs targetVerticalEccDegs], ...
            'rodIntrusionAdjustedConeAperture', true, ...
            'customDegsToMMsConversionFunction', @RGCmodels.Watson.convert.rhoDegsToMMs, ...
            'customMMsToDegsConversionFunction', @RGCmodels.Watson.convert.rhoMMsToDegs);


    % Compute cone indices and weights for the center & surround
    % mechanism based on the retinalRFparamsStruct and the targetConeMosaic
    pooledConeIndicesAndWeightsStruct = weightsComputeFunctionHandle(...
            retinalRFparamsStruct, ...
            targetVisualRFDoGparams.conesNumPooledByTheRFcenter, ...
            [], [], targetConeMosaic);


    % Generate optics for the targetConeMosaic
    switch (opticsParams.ZernikeDataBase)
            % Artal
            case RetinaToVisualFieldTransformer.Artal
                subtractCentralRefraction = ArtalOptics.constants.subjectRequiresCentralRefractionCorrection(...
                    targetConeMosaic.whichEye, opticsParams.subjID);
                % Polans
            case RetinaToVisualFieldTransformer.Polans
                subtractCentralRefraction = PolansOptics.constants.subjectRequiresCentralRefractionCorrection(...
                    targetConeMosaic.whichEye, opticsParams.subjID);
    end

    [oiEnsemble, psfEnsemble] = targetConeMosaic.oiEnsembleGenerate(...
            targetConeMosaic.eccentricityDegs, ...
            'zernikeDataBase', opticsParams.ZernikeDataBase, ...
            'subjectID', opticsParams.subjID, ...
            'pupilDiameterMM', opticsParams.pupilDiameterMM, ...
            'zeroCenterPSF', true, ...
            'subtractCentralRefraction', subtractCentralRefraction, ...
            'wavefrontSpatialSamples', opticsParams.wavefrontSpatialSamples, ...
            'warningInsteadOfErrorForBadZernikeCoeffs', true);

    theOI = oiEnsemble{1};
    thePSFData = psfEnsemble{1};

    % Only keep the maxSpatialSupportDegs portion of the PSF
    idx = find(abs(thePSFData.supportX) < maxSpatialSupportDegs*60);
    idy = find(abs(thePSFData.supportY) < maxSpatialSupportDegs*60);
    thePSFData.supportX = thePSFData.supportX(idx);
    thePSFData.supportY = thePSFData.supportY(idy);
    thePSFData.data = thePSFData.data(idy,idx,:);

    % Compute the retinalRF by summing the weighted cone apertures in the
    % center and surround as specified in the computed pooledConeIndicesAndWeightsStruct
    rfSupportX = rfSpatialSupportDegs(:,1);
    rfSupportY = rfSpatialSupportDegs(:,2);
   
    [theRetinalRFcenter, theRetinalRFsurround] = RetinaToVisualFieldTransformer.generateRFsubregionMapsFromPooledCones(...
       rfSupportX,rfSupportY, targetConeMosaic, pooledConeIndicesAndWeightsStruct);

    % And the full cone-pooling based retinal RF
    theRetinalRF = theRetinalRFcenter - theRetinalRFsurround;

    visualizedProfile = 'LSF';
    maxRF  = max([...
        max(abs(targetVisualRF(:)))...
        max(abs(theRetinalRF(:)))...
        ]);
    maxProfile  = max([...
        max(sum(targetVisualRF,1))...
        max(sum(theRetinalRF,1))...
        ]);
    maxPSF = max(thePSFData.data(:));

    % Wavelengths to visualize
    visualizationWavelengths = [450 500 550 600 650];

    for i = 1:numel(visualizationWavelengths)
        visualizationWavelength = visualizationWavelengths(i);

        hFig = figure(1000+i); clf;
        set(hFig, 'Color', [1 1 1], 'Position', [10 10 1200 950]);
        [~,iWave] = min(abs(visualizationWavelength-thePSFData.supportWavelength));

        % Select which PSF to use
        theWavePSF = squeeze(thePSFData.data(:,:,iWave));
        %theWavePSF = theEmployedPSFData.data;

        % Compute the visualRF
        theWaveVisualRF = conv2(theRetinalRF, theWavePSF, 'same');


        ax = subplot(3,2,1);
        plotRF(ax, rfSupportX, rfSupportY, targetVisualRF, [], maxRF, maxProfile, ...
            visualizedProfile, true, maxSpatialSupportDegs, ...
            sprintf('target visual RF'), true, true);

        ax = subplot(3,2,3);
        plotRF(ax, rfSupportX, rfSupportY, theWaveVisualRF, [], maxRF, ...
            maxProfile, visualizedProfile, true, maxSpatialSupportDegs, ...
            sprintf('achieved @%2.0fnm', visualizationWavelength), true, true);

        ax = subplot(3,2,5);
        plotRF(ax, rfSupportX, rfSupportY, theRetinalRF, [], maxRF, maxProfile, ...
            visualizedProfile, true, maxSpatialSupportDegs, ...
            sprintf('retinal RF\n(derived from V-lambda weighted-\nPSF and target visual RF)'), true, true);


        ax = subplot(3,2,[2 4 6]);
        plotProfiles(ax, thePSFData, theWaveVisualRF, targetVisualRF, ...
            theRetinalRF, visualizedProfile, maxSpatialSupportDegs, ...
            sprintf('achieved @%2.0fnm', visualizationWavelength));

    end



    maxRF  = max([...
        max(abs(targetVisualRF(:)))...
        ]);
    maxProfile  = max([...
        max(sum(targetVisualRF,1))...
        ]);

    hFig = figure(100); clf;
    set(hFig, 'Color', [1 1 1], 'Position', [10 10 1700 750]);

    wavelengthsToExamine = [450 475 500 525 550 575 600 625 650];
    subplotPosVectors = NicePlot.getSubPlotPosVectors(...
        'rowsNum', 3, ...
        'colsNum', numel(wavelengthsToExamine), ...
        'heightMargin',  0.03, ...
        'widthMargin',    0.02, ...
        'leftMargin',     0.01, ...
        'rightMargin',    0.00, ...
        'bottomMargin',   0.01, ...
        'topMargin',      0.01);


    for i = 1:numel(wavelengthsToExamine)

        % Convolve the cone-pooling based retinal RF to get the corresponding visual RF
        [~,iWave] = min(abs(wavelengthsToExamine(i)-thePSFData.supportWavelength));
        theWavePSF = squeeze(thePSFData.data(:,:,iWave));
        theWaveVisualRF = conv2(theRetinalRF, theWavePSF, 'same');

        ax = subplot('Position', subplotPosVectors(1, i).v);
        plotRF(ax, rfSupportX, rfSupportY, theWaveVisualRF, targetVisualRF, maxRF, ...
            maxProfile, visualizedProfile, false, maxSpatialSupportDegs, ...
            sprintf('red: visual RF (%2.0fnm)\nblue: target RF', thePSFData.supportWavelength(iWave)), true, false, ...
            [1 0 0], [0 0 1]);

        ax = subplot('Position', subplotPosVectors(2, i).v);
        plotRF(ax, rfSupportX, rfSupportY, targetVisualRF-theWaveVisualRF, [], maxRF, ...
            maxProfile, visualizedProfile, false, maxSpatialSupportDegs, ...
            sprintf('residual RF\ntarget-achieved (%2.0fnm)', thePSFData.supportWavelength(iWave)), true, false, ...
            [0 0 0], [0 0 0]);

        ax = subplot('Position', subplotPosVectors(3, i).v);
        plotPSF(ax, thePSFData, maxPSF, maxSpatialSupportDegs, iWave);
    end


end

function shadedAreaPlot(ax,x,y, baseline, faceColor, edgeColor, faceAlpha, lineWidth, lineStyle)
    x = [x fliplr(x)];
    y = [y y*0+baseline];
    px = reshape(x, [1 numel(x)]);
    py = reshape(y, [1 numel(y)]);
    pz = -10*eps*ones(size(py)); 
    patch(ax,px,py,pz,'FaceColor',faceColor,'EdgeColor', edgeColor, ...
        'FaceAlpha', faceAlpha, 'LineWidth', lineWidth, 'LineStyle', lineStyle);
end


function plotPSF(ax, thePSFData, maxPSF, maxSpatialSupportDegs, iWave)
    psfZLevels = 0.05:0.1:0.95;
    theWavePSF = thePSFData.data(:,:,iWave);
    contourf(ax,thePSFData.supportX/60, thePSFData.supportY/60, theWavePSF/max(theWavePSF(:)), psfZLevels);
    hold on;
    midRow = (size(thePSFData.data,1)-1)/2+1;
    plot(ax, thePSFData.supportX/60, -maxSpatialSupportDegs*0.75 + 1.7*theWavePSF(midRow,:)/maxPSF*maxSpatialSupportDegs, 'r-', 'LineWidth', 1.5);
    axis(ax,'image'); axis 'xy';


    if (maxSpatialSupportDegs < 0.2)
        tickSeparationDegs = 0.05;
    elseif (maxSpatialSupportDegs < 0.4)
        tickSeparationDegs = 0.1;
    elseif (maxSpatialSupportDegs < 0.6)
        tickSeparationDegs = 0.15;
    else
        tickSeparationDegs = 0.2;
    end

    set(ax, 'XLim', maxSpatialSupportDegs*[-1 1], 'YLim', maxSpatialSupportDegs*[-1 1], ...
        'XTick', -5:tickSeparationDegs:5, 'YTick', -5:tickSeparationDegs:5, 'CLim', [0 1], 'FontSize', 14);
    set(ax, 'XTickLabel', {}, 'YTickLabel', {});
    grid(ax, 'on');
    xlabel(ax,sprintf('%2.2fdegs', 2*maxSpatialSupportDegs));
    title(ax, sprintf('PSF (%2.0fnm)', thePSFData.supportWavelength(iWave)), 'FontWeight','normal');
    colormap(ax,brewermap(1024, 'greys'));
end

function plotRF(ax, rfSupportX, rfSupportY, RF, targetRF, maxRF, maxProfile, visualizedProfile, ...
    renderContourPlot, maxSpatialSupportDegs, titleString, noTickLabels, showColorBar, ...
    achievedRFprofileColor, targetRFprofileColor)
   
    if (renderContourPlot)  
        rfZLevels = -0.91:0.02:0.91;
        contourf(ax,rfSupportX, rfSupportY, RF/maxRF, rfZLevels, 'LineColor', 'none');
    else
        imagesc(ax, rfSupportX, rfSupportY, RF/maxRF);
        hold on;
        switch visualizedProfile
            case 'midRow'
                midRow = (size(RF,1)-1)/2+1;
                theProfile = RF(midRow,:)/maxProfile;
                if (~isempty(targetRF))
                    theTargetProfile = targetRF(midRow,:)/maxProfile;
                end
            case 'LSF'
                theProfile = sum(RF,1)/maxProfile;
                if (~isempty(targetRF))
                    theTargetProfile = sum(targetRF,1)/maxProfile;
                end
        end

    
        plot(ax, rfSupportX, -maxSpatialSupportDegs*0.5 + 1.0*theProfile*maxSpatialSupportDegs, '-', 'Color', achievedRFprofileColor, 'LineWidth', 1.5);
        if (~isempty(targetRF))
            plot(ax, rfSupportX, -maxSpatialSupportDegs*0.5 + 1.0*theTargetProfile*maxSpatialSupportDegs, '-', 'Color', targetRFprofileColor, 'LineWidth', 1.0);
        end
        plot(ax, rfSupportX, -maxSpatialSupportDegs*0.5 + rfSupportX*0, 'k-');
    
    end

    if (maxSpatialSupportDegs < 0.2)
        tickSeparationDegs = 0.05;
    elseif (maxSpatialSupportDegs < 0.4)
        tickSeparationDegs = 0.1;
    elseif (maxSpatialSupportDegs < 0.6)
        tickSeparationDegs = 0.15;
    else
        tickSeparationDegs = 0.2;
    end

    axis(ax,'image'); axis 'xy';
    set(ax, 'XLim', maxSpatialSupportDegs*[-1 1], 'YLim', maxSpatialSupportDegs*[-1 1], ...
        'XTick', -5:tickSeparationDegs:5, 'YTick', -5:tickSeparationDegs:5, 'CLim', [-1 1], 'FontSize', 14);

    if (noTickLabels)
        set(ax, 'XTickLabel', {}, 'YTickLabel', {});
    end
    xtickangle(ax, 90);

    if (showColorBar)
        colorbar(ax, 'EastOutside');
        ylabel(ax,sprintf('space (%2.2fdegs)', 2*maxSpatialSupportDegs));
    end

    grid(ax, 'on');
    colormap(ax,brewermap(1024, '*RdBu'));
    xlabel(ax,sprintf('space (%2.2fdegs)', 2*maxSpatialSupportDegs));
    title(ax, titleString, 'FontWeight','normal');
end


function plotProfiles(ax, rfSupportX, achievedRF, targetRF, retinalRF, visualizedProfile, maxSpatialSupportDegs, achievedLabel)
    
    switch visualizedProfile
        case 'midRow'
            midRow = (size(achievedRF,1)-1)/2+1;
            theAchievedRFProfile = achievedRF(midRow,:);
            theTargetProfile = targetRF(midRow,:);
            theRetinalProfile = retinalRF(midRow,:);
        case 'LSF'
            theAchievedRFProfile = sum(achievedRF,1);
            theTargetProfile = sum(targetRF,1);
            theRetinalProfile = sum(retinalRF,1);
    end
    
    maxProfile = max([max(abs(theAchievedRFProfile)) max(abs(theTargetProfile)) max(abs(theRetinalProfile))]);
    theAchievedRFProfile = theAchievedRFProfile / maxProfile;
    theTargetProfile = theTargetProfile / maxProfile;
    theRetinalProfile = theRetinalProfile / maxProfile;

    faceColor = [0.5 1 0.9];
    edgeColor = faceColor*0.7;
    faceAlpha = 0.8;
    lineWidth = 1.5;
    lineStyle = '-';
    shadedAreaPlot(ax,rfSupportX, theTargetProfile, 0, faceColor, edgeColor, faceAlpha, lineWidth, lineStyle)
    hold(ax, 'on');
    plot(ax, rfSupportX/60, theAchievedRFProfile, 'r-', 'LineWidth', 1.5);
   
    plot(ax, rfSupportX, theTargetProfile-theAchievedRFProfile, 'k--', 'LineWidth', 1.0);
    plot(ax, rfSupportX0, theRetinalProfile, 'k-', 'LineWidth', 1.5);

    if (maxSpatialSupportDegs < 0.2)
        tickSeparationDegs = 0.05;
    elseif (maxSpatialSupportDegs < 0.4)
        tickSeparationDegs = 0.1;
    elseif (maxSpatialSupportDegs < 0.6)
        tickSeparationDegs = 0.15;
    else
        tickSeparationDegs = 0.2;
    end

    set(ax, 'XLim', [min(thePSFData.supportX) max(thePSFData.supportX)]/60, 'YLim', [-0.2 1], ...
            'XTick', -5:tickSeparationDegs:5, 'YTick', -0.6:0.1:1, 'FontSize', 14);
    grid(ax, 'on');
    legend({'target', achievedLabel, sprintf('target - %s', achievedLabel), 'retinal RF'}, 'Location', 'NorthOutside', 'numColumns', 2);
    xlabel(ax,'degrees');
    xtickangle(ax, 90)
end
