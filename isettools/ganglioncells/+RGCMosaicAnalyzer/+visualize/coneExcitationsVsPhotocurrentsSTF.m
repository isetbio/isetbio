%
% RGCMosaicAnalyzer.visualize.coneExcitationsVsPhotocurrentsSTF()
%

function hFig = coneExcitationsVsPhotocurrentsSTF(...
    theConeModulationsBasedSTFamplitudeSpectra, ...
    thePhotocurrentsBasedSTFamplitudeSpectra, ...
    theConeModulationsBasedSTFphaseSpectra, ...
    thePhotocurrentsBasedSTFphaseSpectra, ...
    theConeModulationsBasedResponses, ...
    thePhotocurrentsBasedResponses, ...
    theConeModulationsBasedResponseTemporalSupportSeconds, ...
    thePhotocurrentsBasedResponseTemporalSupportSeconds, ...
    mRGCsOperateOnBackgroundAdaptedPhotocurrents, ...
    stimParams, theRGCindex, theMRGCmosaic, ...
    maxConeModulationResponses, ...
    maxPhotocurrentResponses, ...
    theConeModulationsBasedBPIs, ...
    thePhotocurrentsBasedBPIs, ...
    varargin)


    p = inputParser;
    p.addParameter('debugPhotocurrentPhaseAlignment', false, @islogical);
    p.addParameter('withSuperimposedPSF', [], @(x)(isempty(x)||isstruct(x)));
    p.addParameter('spatialSupportTickSeparationArcMin', 7, @isscalar);
    p.addParameter('extraPlotType', 'differentialPhaseSpectraAndBPIs', @(x)(ismember(x, {'differentialPhaseSpectraAndBPIs', 'individualPhaseSpectra', '2DSTFs'})));
    p.addParameter('exportPDFdirectory', '', @(x)(isempty(x)||ischar(x)));
    p.addParameter('videoOBJ', []);
    % Execute the parser
    p.parse(varargin{:});

    debugPhotocurrentPhaseAlignment = p.Results.debugPhotocurrentPhaseAlignment;
    theSuperimposedPSF = p.Results.withSuperimposedPSF;
    spatialSupportTickSeparationArcMin = p.Results.spatialSupportTickSeparationArcMin;
    extraPlotType = p.Results.extraPlotType;
    exportPDFdirectory = p.Results.exportPDFdirectory;
    videoOBJ = p.Results.videoOBJ;


    coneResponseExtraSamplesNum = 1;
    photocurrentResponseExtraSamplesNum = 0;

    [theTimeAlignedPhotocurrentsResponses, theTimeAlignedConeModulationsResponses, theTimeAlignedTemporalSupportSeconds] = ...
        timeAlignConeModulationAndPhotocurrentResponses(...
            theConeModulationsBasedResponses, thePhotocurrentsBasedResponses, ...
            theConeModulationsBasedResponseTemporalSupportSeconds, ...
            thePhotocurrentsBasedResponseTemporalSupportSeconds, ...
            theConeModulationsBasedSTFphaseSpectra, thePhotocurrentsBasedSTFphaseSpectra, ...
            stimParams, coneResponseExtraSamplesNum, photocurrentResponseExtraSamplesNum, debugPhotocurrentPhaseAlignment);


    % Ticks, ticklabels, and ranges
    [timeTicks, coneModulationsResponseTicks, photocurrentsResponseTicks, ...
     coneModulationsResponseTickLabels, photocurrentsResponseTickLabels, ...
     coneModulationsResponseTicksDouble, photocurrentsResponseTicksDouble, ...
     coneModulationsResponseTickLabelsDouble, photocurrentsResponseTickLabelsDouble, ...
     coneModulationsResponseTicksLog, photocurrentsResponseTicksLog, ...
     coneModulationsResponseTickLogLabels, photocurrentsResponseTickLogLabels, ...
     coneModulationsResponseRange, photocurrentsReponseRange, ...
     coneModulationsResponseRangeLog, photocurrentsReponseRangeLog] = generateTicksAndTickLabels(....
                theConeModulationsBasedResponseTemporalSupportSeconds, ...
                maxConeModulationResponses, maxPhotocurrentResponses, ...
                mRGCsOperateOnBackgroundAdaptedPhotocurrents);

    axesFontSize = 20;
    titleFontSize = 16;


    subplotPosVectors = NicePlot.getSubPlotPosVectors(...
       'rowsNum', 2, ...
       'colsNum', 5, ...
       'heightMargin',  0.11, ...
       'widthMargin',    0.06, ...
       'leftMargin',     0.04, ...
       'rightMargin',    0.00, ...
       'bottomMargin',   0.07, ...
       'topMargin',      0.05);


    hFig = figure(1); clf;
    set(hFig, 'Position', [10 10 2820 1350], 'Color', [1 1 1]);

    ff = PublicationReadyPlotLib.figureComponents('1x1 standard figure',...
        'darkScheme', true);

    if (~isempty(exportPDFdirectory))
        thePDFFilename = sprintf('RGC_%d_nominalC_%2.0f%%_%2.0fCDM2_%2.1fHz.pdf', ...
            theRGCindex, stimParams.contrast*100, stimParams.backgroundLuminanceCdM2, stimParams.temporalFrequencyHz);

        % Generate figure dir if it does not exist
        theScriptName = strrep(exportPDFdirectory, ISETBioPaperAndGrantCodeRootDirectory, '');
        theScriptName = strrep(theScriptName, '/','');
        theFiguresDir = ISETBioPaperAndGrantCodeFigureDirForScript(theScriptName);
    end


    % The cone pooling map & its horizontal line weighting functions
    %axConePoolingMap = subplot('Position', subplotPosVectors(1,1).v);
    axConePoolingMap = subplot('Position', subplotPosVectors(1,5).v);
    %axConePoolingLineWeightingFunctions = subplot('Position', subplotPosVectors(2,1).v);
    axConePoolingLineWeightingFunctions = subplot('Position', subplotPosVectors(2,5).v);
    renderConePoolingMapWithPSF(axConePoolingMap, axConePoolingLineWeightingFunctions, theMRGCmosaic, theRGCindex, ...
        theSuperimposedPSF, spatialSupportTickSeparationArcMin, axesFontSize, titleFontSize);


    % Cone modulation based mRGC responses (zero phase)
    %ax = subplot('Position', subplotPosVectors(1,2).v);
    ax = subplot('Position', subplotPosVectors(1,1).v);

    renderPhaseAlignedResponses(ax, 'stairs plot', stimParams, ...
        theConeModulationsBasedResponseTemporalSupportSeconds, ...
        theConeModulationsBasedResponses, ...
        theConeModulationsBasedSTFphaseSpectra, ...
        coneResponseExtraSamplesNum, coneModulationsResponseRange, [0 theConeModulationsBasedResponseTemporalSupportSeconds(end)], ...
        timeTicks, coneModulationsResponseTicks, coneModulationsResponseTickLabels, ...
        sprintf('cone modulations-based\nmRGC response'), ...
        sprintf('%2.2fHz, %2.0fcd/m2, %2.0f%%', stimParams.temporalFrequencyHz, stimParams.backgroundLuminanceCdM2, 100*stimParams.contrast'), ...
        axesFontSize, titleFontSize, ff);


    % Photocurrent based mRGC responses
    % Aligned with respect to the cone modulation response phase
    % so as to reveal the temporal delay of the photocurrent response
    % with respect to the coneModulation response
    %ax = subplot('Position', subplotPosVectors(2,2).v);
    ax = subplot('Position', subplotPosVectors(1,4).v);
    renderPhaseAlignedResponses(ax, 'line plot', stimParams, ...
        thePhotocurrentsBasedResponseTemporalSupportSeconds, ...
        thePhotocurrentsBasedResponses, ...
        theConeModulationsBasedSTFphaseSpectra, ...
        photocurrentResponseExtraSamplesNum, photocurrentsReponseRange, [0 theConeModulationsBasedResponseTemporalSupportSeconds(end)], ...
        timeTicks, photocurrentsResponseTicks, photocurrentsResponseTickLabels, ...
        'mRGC response (pAmps)', ...
        sprintf('photocurrents-based mRGC responses\n(phase relative to cone modulations response)'), ...
        axesFontSize, titleFontSize, ff);


    % Plotocurrent - based mRGC response time series (zero phase)
    %ax = subplot('Position', subplotPosVectors(1,3).v);
    ax = subplot('Position', subplotPosVectors(1,2).v);

    hold (ax, 'on');
    for iORI = 1:numel(stimParams.orientationDegs)
        if (iORI == 1)
            sfColors = brewermap(5+numel(stimParams.spatialFrequencyCPD), 'reds');
        else
            sfColors = brewermap(5+numel(stimParams.spatialFrequencyCPD), 'blues');
        end
        for iSF = 1:numel(stimParams.spatialFrequencyCPD)
            plot(ax, theTimeAlignedTemporalSupportSeconds*1e3, squeeze(theTimeAlignedPhotocurrentsResponses(iORI,iSF,:)), ...
                '-', 'Color', sfColors(5+iSF,:), 'LineWidth', 1.5);
        end
    end
    grid(ax, 'on')
    axis(ax, 'square')
    timeLims = [0 theConeModulationsBasedResponseTemporalSupportSeconds(end)]*1e3;
    set(ax, 'FontSize', 20, 'Ylim', photocurrentsReponseRange, 'XLim', timeLims);
    xtickangle(ax, 0)
    set(ax, 'XTick', timeTicks*1e3, 'YTick', photocurrentsResponseTicks);
    xlabel(ax,'time (msec)',  'FontAngle', 'italic');
    ylabel(ax, sprintf('photocurrents-based\nmRGC response (pAmps)'),  'FontAngle', 'italic');
    title(ax, sprintf('%2.2fHz, %2.0fcd/m2, %2.0f%%', stimParams.temporalFrequencyHz, stimParams.backgroundLuminanceCdM2, 100*stimParams.contrast'), 'FontWeight', 'normal', 'FontSize', 16);

    PublicationReadyPlotLib.offsetAxes(ax,ff, timeLims, photocurrentsReponseRange);
    PublicationReadyPlotLib.applyFormat(ax,ff);


    % The relationship between the time-aligned cone excitations and time-aligned photocurrents
    %ax = subplot('Position', subplotPosVectors(2,3).v);
    ax = subplot('Position', subplotPosVectors(1,3).v);
    renderConeModulationPhotocurrentResponseFunction(ax, stimParams, ...
        theTimeAlignedConeModulationsResponses, theTimeAlignedPhotocurrentsResponses, ...
        coneModulationsResponseRange, photocurrentsReponseRange, ...
        coneModulationsResponseTicks, photocurrentsResponseTicks, ...
        coneModulationsResponseTickLabels, photocurrentsResponseTickLabels, ...
        sprintf('photocurrent non-linearity'), ...
        axesFontSize, titleFontSize, ff);


    % Cone-based STF amplitude spectra for the examined orientations
    %ax = subplot('Position', subplotPosVectors(1,4).v);
    ax = subplot('Position', subplotPosVectors(2,1).v);

    renderSTFamplitudeSpectra(ax, stimParams, ...
        theConeModulationsBasedSTFamplitudeSpectra, ...
        coneModulationsResponseRange,  coneModulationsResponseTicksDouble, coneModulationsResponseTickLabelsDouble, ...
        theConeModulationsBasedBPIs, ...
        '|STF| (modulation)', ...
        sprintf('cone modulations - based'), ...
        axesFontSize, titleFontSize, ff);


    % Photocurrents-based STF amplitude spectra for the examined orientations
    %ax = subplot('Position', subplotPosVectors(1,5).v);
    ax = subplot('Position', subplotPosVectors(2,2).v);

    renderSTFamplitudeSpectra(ax, stimParams, ...
        thePhotocurrentsBasedSTFamplitudeSpectra, ...
        photocurrentsReponseRange, photocurrentsResponseTicksDouble, photocurrentsResponseTickLabelsDouble, ...
        thePhotocurrentsBasedBPIs, ...
        '|STF| (pAmps)', ...
        sprintf('photocurrents - based'), ...
        axesFontSize, titleFontSize, ff);



    switch (extraPlotType)
        case 'differentialPhaseSpectraAndBPIs'
            
            centerPhaseDegs = 90;
            phaseRangeDegs = 120;

            ax = subplot('Position', subplotPosVectors(1,4).v);
            renderDifferentialSTFphaseSpectra(ax, stimParams, ...
                theConeModulationsBasedSTFphaseSpectra, thePhotocurrentsBasedSTFphaseSpectra, ...
                centerPhaseDegs, phaseRangeDegs, ...
                sprintf('differential STF phase spectra\n(photocurrents - cone modulations)'), ...
                axesFontSize, titleFontSize);

            % The BPI's for all cells up to this one
            %ax = subplot('Position', subplotPosVectors(2,5).v);
            ax = subplot('Position', subplotPosVectors(2,3).v);
            renderBPIcorrespondence(ax, stimParams, ...
                theConeModulationsBasedBPIs, thePhotocurrentsBasedBPIs, ...
                axesFontSize, titleFontSize, ff);
            PublicationReadyPlotLib.applyFormat(ax,ff);
           

        case 'individualPhaseSpectra'
            % Cone-based STF phase spectra for the examined orientations
            ax = subplot('Position', subplotPosVectors(1,4).v);
            renderSTFphaseSpectra(ax, stimParams, ...
                theConeModulationsBasedSTFphaseSpectra, ...
                'cone modulations - based STF phase', axesFontSize, titleFontSize);
        

            % Photocurrents-based STF phase spectra for the examined orientations
            ax = subplot('Position', subplotPosVectors(2,5).v);
            renderSTFphaseSpectra(ax, stimParams, ...
                thePhotocurrentsBasedSTFphaseSpectra, ...
                'photocurrents - based STF phase', axesFontSize, titleFontSize);

        case '2DSTFs'

            [coneModulationsSTF2D, coneModulationsSTFMatrix] = generate2DSTF(stimParams, theConeModulationsBasedSTFamplitudeSpectra);
            [photocurrentsSTF2D, photocurrentsSTFMatrix] = generate2DSTF(stimParams, thePhotocurrentsBasedSTFamplitudeSpectra);

            % 2D STF for cone modulations
            ax = subplot('Position', subplotPosVectors(1,4).v);
            render2DSTFamplitudeSpectra(ax, ...
                coneModulationsSTFMatrix, coneModulationsSTF2D, ...
                sfTicks, sfTickLabels, ...
                'cone modulations-based 2D STF', axisFontSize, titleFontSize);

            % 2D STF for photocurrents
            ax = subplot('Position', subplotPosVectors(2,5).v);
            render2DSTFamplitudeSpectra(ax, ...
                photocurrentsSTFMatrix, photocurrentsSTF2D, ...
                sfTicks, sfTickLabels, ...
                'photocurrents-based 2D STF', axisFontSize, titleFontSize);

        otherwise
            error('Unknown extraPlotType: ''%s''.', extraPlotType);

    end


    drawnow;

    if (~isempty(videoOBJ))
        videoOBJ.writeVideo(getframe(hFig));
    end

    if (~isempty(exportPDFdirectory))
        NicePlot.exportFigToPDF(fullfile(theFiguresDir,thePDFFilename), hFig, 300);
    end
end


function  [theSTF2D, theSTFMatrix] = generate2DSTF(stimParams, theSTFamplitudeSpectra)

    x = [];
    y = [];
    theSTFMatrix = [];

    for iSF = 1:numel(stimParams.spatialFrequencyCPD)
        for iORI = 1:numel(stimParams.orientationDegs)
            x(numel(x)+1) = iSF * cosd(stimParams.orientationDegs(iORI));
            y(numel(y)+1) = iSF * sind(stimParams.orientationDegs(iORI));
            theSTFMatrix(numel(theSTFMatrix)+1) = theSTFamplitudeSpectra(iORI,iSF);

            x(numel(x)+1) = -x(numel(x));
            y(numel(y)+1) = -y(numel(y));
            theSTFMatrix(numel(theSTFMatrix)+1) = theSTFamplitudeSpectra(iORI,iSF);
        end
    end

    interpolationMethod = 'natural';
    extrapolationMethod = 'none';
    F = scatteredInterpolant(x(:), y(:), theSTFMatrix(:), interpolationMethod, extrapolationMethod);

    xx = generateSFticks(stimParams.spatialFrequencyCPD);

    yy = xx;
    [X,Y] = meshgrid(xx,yy);
    theSTF2D = F(X(:),Y(:));
    theSTF2D = reshape(theSTF2D , numel(xx), numel(yy));
end

function render2DSTFamplitudeSpectra(ax, ...
                theSTFMatrix, theSTF2D, ...
                sfTicks, sfTickLabels, ...
                plotTitle, axisFontSize, titleFontSize)

    %imagesc(ax, xx,yy,theSTF2D);

    zLevels = (0:0.1:1.0)*max(theSTFMatrix(:));
    contourf(ax, xx, yy, theSTF2D, zLevels);

    set(ax, 'XTick', [sfTicks(1) 0 sfTicks(end)], 'XTickLabel', {sfTickLabels{1}, '0', sfTickLabels{end}}, ...
        'YTick', sfTicks, 'YTickLabel', sfTickLabels);
    colormap(ax, brewermap(1024, '*greys'));
    axis(ax, 'image');
    set(ax, 'CLim', [0 max(theSTFMatrix(:))], 'Color', [0 0 0]);
    set(ax, 'FontSize', titleFontSize)
    xtickangle(ax, 0)
    title(ax, plotTitle, 'FontWeight', 'normal', 'FontSize', axisFontSize);
    xlabel(ax,'spatial frequency (c/deg)',  'FontAngle', 'italic');
    ylabel(ax,'spatial frequency (c/deg)',  'FontAngle', 'italic');
end


function renderBPIcorrespondence(ax, stimParams, ...
                theConeModulationsBasedBPIs, thePhotocurrentsBasedBPIs, ...
                axesFontSize, titleFontSize, ff)
            
    hold(ax,'on')
    theLegends = cell(1,numel(stimParams.orientationDegs));

    

    lineColors1 = brewermap(5, 'reds');
    lineColors2 = brewermap(5, 'blues');


    oriColors = [lineColors1(5,:); lineColors2(5,:)];

    pHandles = [];

    % All previous cell BPIs (small disks)
    for iCell = 1:size(theConeModulationsBasedBPIs,1)
          plot(ax,theConeModulationsBasedBPIs(iCell, :), thePhotocurrentsBasedBPIs(iCell, :), 'k-', 'LineWidth', 3.0);
          plot(ax,theConeModulationsBasedBPIs(iCell, :), thePhotocurrentsBasedBPIs(iCell, :), '-','Color', [0.5 0.5 0.5], 'LineWidth', 1.0);
    end

    for iORI = 1:numel(stimParams.orientationDegs) 
        pHandles(iORI) = scatter(ax,theConeModulationsBasedBPIs(:, iORI), thePhotocurrentsBasedBPIs(:, iORI), 16^2, ...
            'MarkerFaceColor', squeeze(oriColors(iORI,:)), 'MarkerFaceAlpha', 1.0, ...
            'MarkerEdgeColor', squeeze(oriColors(iORI,:))*0.5, 'LineWidth', 1.5, 'MarkerEdgeAlpha', 1.0);
        theLegends{iORI} = sprintf('%d degs', stimParams.orientationDegs(iORI));
    end

    

    outlineCurrentCellInGreen = false;
    if (outlineCurrentCellInGreen)
        % Black line joining the BPIs for all orientations for the current cell
        plot(ax,theConeModulationsBasedBPIs(end, :), thePhotocurrentsBasedBPIs(end, :), 'k-', 'LineWidth', 3.0);
        plot(ax,theConeModulationsBasedBPIs(end, :), thePhotocurrentsBasedBPIs(end, :), '-','Color', 'g', 'LineWidth', 1.0);
    
        % Current cell BPIs (large cyan square)
        for iORI = 1:numel(stimParams.orientationDegs)
            scatter(ax,theConeModulationsBasedBPIs(end, iORI), thePhotocurrentsBasedBPIs(end, iORI), 20^2, ...
                'Marker', 's', 'MarkerFaceColor', [0.0 1 0.0] , 'MarkerFaceAlpha', 1.0, ...
                'MarkerEdgeColor', [0 0 0], 'MarkerEdgeAlpha', 1.0, 'LineWidth', 2.0);
        end
    end

    plot([0 1], [0 1], 'w-', 'LineWidth', 1.5);

    axis(ax, 'square');
    set(ax, 'XTick', 0:0.2:1, 'XTickLabel', {'0.0', '0.2', '0.4', '0.6', '0.8', 'LP'});
    set(ax, 'YTick', 0:0.2:1, 'YTickLabel', {'0.0', '0.2', '0.4', '0.6', '0.8', 'LP'});
    set(ax, 'XLim', [0 1], 'YLim', [0 1]);
    grid(ax, 'on');
    legend(ax, pHandles, theLegends, 'Location', 'SouthEast');
    set(ax, 'FontSize', axesFontSize)
    xtickangle(ax, 0)
    xlabel(ax,'cone modulations-based', 'FontAngle', 'italic');
    ylabel(ax,'photocurrents-based',  'FontAngle', 'italic');
    title(ax, sprintf('STF bandpass index\n|STF(0)| / |STF(peak sf)|'), 'FontWeight', 'normal', 'FontSize', titleFontSize);

    PublicationReadyPlotLib.offsetAxes(ax,ff, [0 1], [0 1]);
    PublicationReadyPlotLib.applyFormat(ax,ff);
end


function renderSTFphaseSpectra(ax, stimParams, theSTFphaseSpectra, ...
                theYaxisLabel, axesFontSize, titleFontSize)
                

    hold(ax,'on')
    theLegends = cell(1,numel(stimParams.orientationDegs));

    oriColors = brewermap(numel(stimParams.orientationDegs), 'RdBu');

    % Tolerance for a phase jump between consecutive angles that when we exeed,
    % we shift the angles by adding multiples of 360 degs
    toleranceDegs = 180;
    toleranceRadians = toleranceDegs/180*pi;

    for iORI = 1:numel(stimParams.orientationDegs)
        theWrappedPhaseDegs = squeeze(theSTFphaseSpectra(iORI,:));
        theUnwrappedPhaseDegs = unwrap(theWrappedPhaseDegs/180*pi, toleranceRadians)/pi*180;
        plot(ax, stimParams.spatialFrequencyCPD, theUnwrappedPhaseDegs, ...
            'o-', 'Color', squeeze(oriColors(iORI,:)), 'LineWidth', 1.5, ...
            'MarkerSize', 12, 'MarkerFaceColor', squeeze(oriColors(iORI,:)), 'MarkerEdgeColor', 0.5*squeeze(oriColors(iORI,:)));
        theLegends{iORI} = sprintf('%d degs', stimParams.orientationDegs(iORI));
    end
    axis(ax, 'square')
    yTicks = -360:90:360;
    set(ax, 'XScale', 'log', 'XTick', [0.01 0.03 0.1 0.3 1 3 10 30]);
    set(ax, 'YTick', yTicks);
    set(ax, 'XLim', [0.05 60], 'YLim', [-360 360]);
    grid(ax, 'on');
    legend(ax, theLegends, 'Location', 'SouthWest');
    set(ax, 'FontSize', axesFontSize)
    xtickangle(ax, 0)
    xlabel(ax,'spatial frequency (c/deg)',  'FontAngle', 'italic');
    ylabel(ax, theYaxisLabel,  'FontAngle', 'italic');
end

function renderDifferentialSTFphaseSpectra(ax, stimParams, ...
    theConeModulationsBasedSTFphaseSpectra, thePhotocurrentsBasedSTFphaseSpectra, ...
    centerPhaseDegs, phaseRangeDegs, ...
    plotTitle, axesFontSize, titleFontSize)

    hold(ax,'on')
    theLegends = cell(1,numel(stimParams.orientationDegs));
    
    oriColors = brewermap(numel(stimParams.orientationDegs), 'RdBu');

    % Tolerance for a phase jump between consecutive angles that when we exeed,
    % we shift the angles by adding multiples of 360 degs
    toleranceDegs = 180;
    toleranceRadians = toleranceDegs/180*pi;

    for iORI = 1:numel(stimParams.orientationDegs)
        if (iORI == 1)
            lineColors = brewermap(5+numel(stimParams.spatialFrequencyCPD), 'reds');
        else
            lineColors = brewermap(5+numel(stimParams.spatialFrequencyCPD), 'blues');
        end
    
        theConeModulationBasedSTFphaseDegs = squeeze(theConeModulationsBasedSTFphaseSpectra(iORI,:));
        thePhotocurrentBasedSTFphaseDegs = squeeze(thePhotocurrentsBasedSTFphaseSpectra(iORI,:));
    
        % To radians
        theConeModulationBasedSTFphaseRadians = theConeModulationBasedSTFphaseDegs/180*pi;
        thePhotocurrentBasedSTFphaseRadians = thePhotocurrentBasedSTFphaseDegs/180*pi;
    
        % To get the phase difference, we divide the Euler form of the normalized STFs
        % This is because ratios of complex numbers in Euler form (r * exp(i * theta)
        % involves dividing their magnitudes (r) and subtracting their angles (phase)
    
        EulerFormNormalizedPhotocurrentSTFs = exp(1j*thePhotocurrentBasedSTFphaseRadians);
        EulerFormNormalizedConeModulationSTFs = exp(1j*theConeModulationBasedSTFphaseRadians);
        theEulerRatio = EulerFormNormalizedPhotocurrentSTFs ./ EulerFormNormalizedConeModulationSTFs;
        thePhaseDifferenceRadians = angle(theEulerRatio);
    
        % Unwrap and back to degrees
        theUnwrappedPhaseDifferenceDegs = unwrap(thePhaseDifferenceRadians, toleranceRadians)/pi*180;
    
        plot(ax, stimParams.spatialFrequencyCPD, theUnwrappedPhaseDifferenceDegs, ...
            '-', 'Color', squeeze(lineColors(5+7,:)), 'LineWidth', 1.5);
    
        for iSF = 1:numel(stimParams.spatialFrequencyCPD)
            p = plot(ax,stimParams.spatialFrequencyCPD(iSF), squeeze(theUnwrappedPhaseDifferenceDegs(iSF)), ...
                'o-', 'Color', squeeze(lineColors(5+iSF,:)), 'LineWidth', 1.5, ...
                'MarkerSize', 12, 'MarkerFaceColor', squeeze(lineColors(5+iSF,:)), 'MarkerEdgeColor', 0.5*squeeze(lineColors(5+iSF,:)));
            if (iSF == 7)
                thePlotHandles(iORI) = p;
                theLegends{iORI} = sprintf('%d degs', stimParams.orientationDegs(iORI));
            end
        end
    end

    axis(ax, 'square')
    yTicks = -360:15:360;
    set(ax, 'XScale', 'log', 'XTick', [0.01 0.03 0.1 0.3 1 3 10 30]);
    set(ax, 'YTick', yTicks);
    set(ax, 'XLim', [0.05 60], 'YLim', centerPhaseDegs + phaseRangeDegs*[-0.5 0.5]);
    grid(ax, 'on');
    legend(ax, thePlotHandles, theLegends, 'Location', 'SouthEAST');
    set(ax, 'FontSize', axesFontSize)
    xtickangle(ax, 0)
    xlabel(ax,'spatial frequency (c/deg)',  'FontAngle', 'italic');
    ylabel(ax,'STF phase difference (degs)',  'FontAngle', 'italic');
    title(ax, plotTitle, 'FontWeight', 'normal', 'FontSize', titleFontSize);

end

function renderSTFamplitudeSpectra(ax, stimParams, ...
        theSTFamplitudeSpectra, ...
        theResponseRange, theResponseTicks, theResponseTickLabels, ...
        theBPIs, ...
        theYaxisLabel, ...
        plotTitle, ...
        axesFontSize, titleFontSize, ff)

    hold(ax,'on')

    thePlotHandles = [];
    for iORI = 1:numel(stimParams.orientationDegs)
        if (iORI == 1)
            lineColors = brewermap(5+numel(stimParams.spatialFrequencyCPD), 'reds');
        else
            lineColors = brewermap(5+numel(stimParams.spatialFrequencyCPD), 'blues');
        end
        plot(ax,stimParams.spatialFrequencyCPD, squeeze(theSTFamplitudeSpectra(iORI,:)), ...
                '-', 'Color', squeeze(lineColors(5+round(0.5*numel(stimParams.spatialFrequencyCPD)),:)), 'LineWidth', 2);

        for iSF = 1:numel(stimParams.spatialFrequencyCPD)
            p = plot(ax,stimParams.spatialFrequencyCPD(iSF), squeeze(theSTFamplitudeSpectra(iORI,iSF)), ...
                'o-', 'Color', squeeze(lineColors(5+iSF,:)), 'LineWidth', 2, ...
                'MarkerSize', 16, 'MarkerFaceColor', squeeze(lineColors(5+iSF,:)), 'MarkerEdgeColor', 0.5*squeeze(lineColors(5+iSF,:)));
            if (iSF == round(0.5*numel(stimParams.spatialFrequencyCPD)))
                thePlotHandles(iORI) = p;
                theLegends{iORI} = sprintf('%d degs (BPI:%2.2f)', stimParams.orientationDegs(iORI), theBPIs(end, iORI));
            end
        end
        
    end
    

    yLims = [0 theResponseRange(2)];
    sfLims = [0.05 60];
    axis(ax, 'square')
    set(ax, 'XScale', 'log', 'XTick', [0.01 0.03 0.1 0.3 1 3 10 30]);
    set(ax, 'YScale', 'linear', 'YTick', theResponseTicks, 'YTickLabel', theResponseTickLabels);
    set(ax, 'XLim', [0.05 60], 'YLim', yLims);
    set(ax, 'YLim', [0 1.1*max(theSTFamplitudeSpectra(:))]);


    grid(ax, 'on');
    legend(ax, thePlotHandles, theLegends, 'Location', 'SouthWest');
    set(ax, 'FontSize', axesFontSize);
    xtickangle(ax, 0)
    xlabel(ax,'spatial frequency (c/deg)',  'FontAngle', 'italic');
    ylabel(ax, theYaxisLabel,  'FontAngle', 'italic');
    title(ax, plotTitle, 'FontWeight', 'normal', 'FontSize', titleFontSize);

    PublicationReadyPlotLib.offsetAxes(ax,ff, sfLims, yLims);
    PublicationReadyPlotLib.applyFormat(ax,ff);
end


function renderPhaseAlignedResponses(ax, thePlotType, stimParams, ...
        theResponseTemporalSupportSeconds, ...
        theResponsesTimeSeries, ...
        theSTFphaseSpectra, ...
        theResponseExtraSamplesNum, theResponseRange, theTemporalSupportRange, ...
        theTimeTicks, theResponseTicks, theResponseTickLabels, ...
        theYaxisLabel, theTitle, axesFontSize, titleFontSize, ff)


    hold (ax, 'on');
    for iORI = 1:numel(stimParams.orientationDegs)
        if (iORI == 1)
            sfColors = brewermap(5+numel(stimParams.spatialFrequencyCPD), 'reds');
        else
            sfColors = brewermap(5+numel(stimParams.spatialFrequencyCPD), 'blues');
        end
        for iSF = 1:numel(stimParams.spatialFrequencyCPD)

            theResponse = squeeze(theResponsesTimeSeries(iORI, iSF, :));
            idx = 1:(numel(theResponse)-theResponseExtraSamplesNum);

            % Phase alignment
            phaseForAlignment = theSTFphaseSpectra(iORI, iSF);
            thePhaseAlignedResponse = phaseAlignResponse(theResponse(idx),...
                phaseForAlignment, ...
                theResponseTemporalSupportSeconds(idx), ...
                1./(stimParams.temporalFrequencyHz), ...
                false);

            if (strcmp(thePlotType, 'stairs plot'))
                stairs(ax, theResponseTemporalSupportSeconds(idx)*1e3, thePhaseAlignedResponse, ...
                    'Color', sfColors(5+iSF,:), 'LineWidth', 1.5);
            else
                plot(ax, theResponseTemporalSupportSeconds(idx)*1e3, thePhaseAlignedResponse, ...
                    '-', 'Color', sfColors(5+iSF,:), 'LineWidth', 1.5);
            end

        end
    end
    grid(ax, 'on')
    axis(ax, 'square')
    set(ax, 'FontSize', axesFontSize, 'Ylim', theResponseRange, 'YTick', theResponseTicks, 'YTickLabel', theResponseTickLabels);
    set(ax, 'XLim', theTemporalSupportRange*1e3);
    set(ax, 'XTick', theTimeTicks*1e3);
    
    xtickangle(ax, 0)
    xlabel(ax,'time (msec)',  'FontAngle', 'italic');
    ylabel(ax, theYaxisLabel,  'FontAngle', 'italic');
    title(ax, theTitle, 'FontWeight', 'normal', 'FontSize', titleFontSize);

    PublicationReadyPlotLib.offsetAxes(ax,ff,theTemporalSupportRange*1e3, theResponseRange);
    PublicationReadyPlotLib.applyFormat(ax,ff);
end


function renderConeModulationPhotocurrentResponseFunction(ax, stimParams, ...
        theTimeAlignedConeModulationsResponses, theTimeAlignedPhotocurrentsResponses, ...
        coneModulationsResponseRange, photocurrentsReponseRange, ...
        coneModulationsResponseTicks, photocurrentsResponseTicks, ...
        coneModulationsResponseTickLabels, photocurrentsResponseTickLabels, ...
        plotTitle, ...
        axesFontSize, titleFontSize, ff)

    lineAlpha = 0.8;
    lineWidth = 1.5;

    for iORI = 1:numel(stimParams.orientationDegs)
        if (iORI == 1)
            lineColors = brewermap(5+numel(stimParams.spatialFrequencyCPD), 'reds');
        else
            lineColors = brewermap(5+numel(stimParams.spatialFrequencyCPD), 'blues');
        end
        for iSF = 1:numel(stimParams.spatialFrequencyCPD)
            theTimeAlignedConeModulationsResponse = theTimeAlignedConeModulationsResponses(iORI,iSF,:);
            theTimeAlignedPhotocurrentsResponse = theTimeAlignedPhotocurrentsResponses(iORI,iSF,:);
            p = patch(ax, [theTimeAlignedConeModulationsResponse(:);NaN],[theTimeAlignedPhotocurrentsResponse(:);NaN],'k');
            set(p, 'edgeColor', lineColors(5+iSF,:), 'edgeAlpha', lineAlpha, 'lineWidth', lineWidth);
        end
    end

    hold(ax, 'on');

    crossHairs = false;
    if (crossHairs)
        plot(ax, [0 0], photocurrentsReponseRange, 'k-', 'LineWidth', 1.0);
        plot(ax, coneModulationsResponseRange, [0 0], 'k-', 'LineWidth', 1.0);
    end

    hold(ax, 'off')
    grid(ax, 'on')
    axis(ax, 'square');
    set(ax, 'FontSize', axesFontSize, ...
        'Ylim', photocurrentsReponseRange, 'XLim', coneModulationsResponseRange, ...
        'XTick', coneModulationsResponseTicks, 'XTickLabel', coneModulationsResponseTickLabels, ...
        'YTick', photocurrentsResponseTicks, 'YTickLabel', photocurrentsResponseTickLabels);
    xtickangle(ax, 0)
    xlabel(ax,sprintf('cone modulations-based\nmRGC response'),  'FontAngle', 'italic');
    ylabel(ax,sprintf('photocurrents-based\nmRGC response (pAmps)'),  'FontAngle', 'italic');
    title(ax, plotTitle, 'FontWeight', 'normal', 'FontSize', titleFontSize);

    
    PublicationReadyPlotLib.offsetAxes(ax,ff,coneModulationsResponseRange, photocurrentsReponseRange);
    PublicationReadyPlotLib.applyFormat(ax,ff);

end


function renderConePoolingMapWithPSF(axConePoolingMap, axConePoolingLineWeightingFunctions, theMRGCmosaic, ...
    theRGCindex, theSuperimposedPSF, spatialSupportTickSeparationArcMin, axesFontSize, titleFontSize)

    surroundConeSelection = 'surround pooling weights > center pooling weights';
    [surroundConePurity, centerConeDominance, centerConeNumerosity, centerConePurity] = ...
			theMRGCmosaic.surroundConePurities(theRGCindex, surroundConeSelection);

    % 10% which is the noise floor in the measurements of cone weights according to RF center overlap according to Greg Field
    % and which is what they used in Fig 4 of their 2010 paper.
    minCenterConeWeightForConePoolingMapVisualization = mRGCMosaic.minRFcenterConeWeightIncludedToMatchFigure4OfFieldEtAl2010;

    % Include surround cones whose pooling weights are >= 0.005 
    % (This is the threshold used in the Field et al 2010 paper -
    % Greg Field - personal communication)
    minSurroundConeWeightForConePoolingMapVisualization = 0.005;

    [scaleBarDegsDefault, ...
     scaleBarMicronsDefault, ...
	 spatialSupportTickSeparationArcMinDefault, ...
	 spatialSupportCenterDegsDefault, ...
     domainVisualizationLimitsDefault, ...
     domainVisualizationTicksDefault, ...
     domainVisualizationLimitsSingleRFDefault, ...
     domainVisualizationTicksSingleRFDefault] = ...
		 	RGCMosaicAnalyzer.visualize.generateLimits(theMRGCmosaic, theMRGCmosaic.rgcRFpositionsDegs(theRGCindex,:), ...
            'spatialSupportTickSeparationArcMin', spatialSupportTickSeparationArcMin);

    scaleBarDegs = scaleBarDegsDefault;

    [~, ~, centerLineWeightingFunctions, surroundLineWeightingFunctions] = ...
        theMRGCmosaic.visualizeCenterSurroundConePoolingMap(theRGCindex, ...
            'axesToRenderIn', axConePoolingMap, ...
            'minConeWeightForVisualizingRFcenterPooling', minCenterConeWeightForConePoolingMapVisualization, ...
            'minConeWeightForVisualizingRFsurroundPooling', minSurroundConeWeightForConePoolingMapVisualization, ...
            'minSurroundConeWeightRelativity', 'center', ...
            'withLineWeightingFunctions', false, ...
            'withSuperimposedPSF', theSuperimposedPSF, ...
            'domainVisualizationLimits', domainVisualizationLimitsSingleRFDefault, ...
            'domainVisualizationTicks', domainVisualizationTicksSingleRFDefault, ...
            'spatialSupportTickSeparationArcMin', spatialSupportTickSeparationArcMin, ...
            'scaleBarDegs', scaleBarDegs, ...
            'doNotLabelScaleBar', true, ...
            'figNo', [], ...
            'figPos', [], ...
            'exportVisualizationPDF',false, ...
            'exportVisualizationPNG',false);
    set(axConePoolingMap, 'FontSize', axesFontSize);
    title(axConePoolingMap, sprintf('RGC #%d\npurities: %2.2f (center), %2.2f (surround)', theRGCindex,  centerConePurity, surroundConePurity), 'FontSize', titleFontSize);
    xtickangle(axConePoolingMap, 0)


    % The horizontal profile of the cone pooling map
    whichMeridian = 'horizontal';		% choose between {'horizontal', 'vertical'}
    spatialSupportCenterDegs = theMRGCmosaic.rgcRFpositionsDegs(theRGCindex,:);
	RGCMosaicConstructor.visualize.centerAndSurroundConePoolingLineWeightingFunctions(...
            '', [], [], ...
		    spatialSupportCenterDegs, spatialSupportTickSeparationArcMin, ...
		    centerLineWeightingFunctions, surroundLineWeightingFunctions, whichMeridian, ...
            'domainVisualizationLimits', domainVisualizationLimitsSingleRFDefault(1:2), ...
            'domainVisualizationTicks', domainVisualizationTicksSingleRFDefault.x, ...
            'yLimsRange', [-0.8 1.01], ...
            'axesToRenderIn', axConePoolingLineWeightingFunctions);
    set(axConePoolingLineWeightingFunctions, 'FontSize', axesFontSize);
    xtickangle(axConePoolingLineWeightingFunctions, 0)
end


function [timeTicks, coneModulationsResponseTicks, photocurrentsResponseTicks, ...
        coneModulationsResponseTickLabels, photocurrentsResponseTickLabels, ...
        coneModulationsResponseTicksDouble, photocurrentsResponseTicksDouble, ...
        coneModulationsResponseTickLabelsDouble, photocurrentsResponseTickLabelsDouble, ...
        coneModulationsResponseTicksLog, photocurrentsResponseTicksLog, ...
        coneModulationsResponseTickLogLabels, photocurrentsResponseTickLogLabels, ... 
        coneModulationsResponseRange, photocurrentsReponseRange, ...
        coneModulationsResponseRangeLog, photocurrentsReponseRangeLog] = generateTicksAndTickLabels(....
                theConeModulationsBasedResponseTemporalSupportSeconds, ...
                maxConeModulationResponses, maxPhotocurrentResponses, ...
                mRGCsOperateOnBackgroundAdaptedPhotocurrents)

    if (theConeModulationsBasedResponseTemporalSupportSeconds(end) <= 0.1)
        timeTicks =  0:0.02:2;
    elseif (theConeModulationsBasedResponseTemporalSupportSeconds(end) <= 0.25)
        timeTicks =  0:0.04:2;
    elseif (theConeModulationsBasedResponseTemporalSupportSeconds(end) <= 0.4)
        timeTicks =  0:0.05:2;
    elseif (theConeModulationsBasedResponseTemporalSupportSeconds(end) <= 0.5)
        timeTicks =  0:0.1:2;
    elseif (theConeModulationsBasedResponseTemporalSupportSeconds(end) <= 1.0)
        timeTicks =  0:0.2:2;
    elseif (theConeModulationsBasedResponseTemporalSupportSeconds(end) <= 3.0)
        timeTicks =  0:0.5:5;
    elseif (theConeModulationsBasedResponseTemporalSupportSeconds(end) <= 5)
        timeTicks =  0:1:5;
    end

    if (maxConeModulationResponses <= 0.2)
        coneModulationsResponseTicks = -0.2:0.1:0.2;
        coneModulationsResponseTicksDouble =  -0.2:0.05:0.2;
        coneModulationsResponseTicksLog = [0.003 0.01 0.03 0.1 0.3];
        coneModulationsResponseTicks = coneModulationsResponseTicks(abs(coneModulationsResponseTicks)<=maxConeModulationResponses);
        coneModulationsResponseTickLabels = sprintf('%+.2f\n', coneModulationsResponseTicks);
        coneModulationsResponseTickLabelsDouble = sprintf('%+.2f\n', coneModulationsResponseTicksDouble);
        coneModulationsResponseTickLogLabels = {'.003', '.01', '.03', '.1', '.3'};
        coneModulationsResponseRange = maxConeModulationResponses*[-1 1];
        coneModulationsResponseRangeLog = [0.001 1.05] * maxConeModulationResponses;
        
    elseif (maxConeModulationResponses <= 0.6)
        coneModulationsResponseTicks = -0.6:0.2:0.6;
        coneModulationsResponseTicksDouble =  -0.6:0.1:0.6;
        coneModulationsResponseTicksLog = [0.003 0.01 0.03 0.1 0.3];
        coneModulationsResponseTicks = coneModulationsResponseTicks(abs(coneModulationsResponseTicks)<=maxConeModulationResponses);
        coneModulationsResponseTickLabels = strrep(sprintf('%+.1f\n', coneModulationsResponseTicks), '0.', '.');
        coneModulationsResponseTickLabelsDouble = strrep(sprintf('%+.1f\n', coneModulationsResponseTicksDouble), '0.', '.');
        coneModulationsResponseTickLogLabels = {'.003', '.01', '.03', '.1', '.3'};
        coneModulationsResponseRange = maxConeModulationResponses*[-1 1];
        coneModulationsResponseRangeLog = [0.001 1.05] * maxConeModulationResponses;

    elseif (maxConeModulationResponses <= 0.8)
        coneModulationsResponseTicks = -0.8:0.4:0.8;
        coneModulationsResponseTicksDouble = -0.8:0.2:0.8;
        coneModulationsResponseTicksLog = [0.003 0.01 0.03 0.1 0.3];
        coneModulationsResponseTicks = coneModulationsResponseTicks(abs(coneModulationsResponseTicks)<=maxConeModulationResponses);
        coneModulationsResponseTickLabels = strrep(strrep(sprintf('%+.1f\n', coneModulationsResponseTicks), '0.', '.'), '1.0', '1');
        coneModulationsResponseTickLabelsDouble = strrep(strrep(sprintf('%+.1f\n', coneModulationsResponseTicksDouble), '0.', '.'), '1.0', '1');
        coneModulationsResponseTickLogLabels = {'.003', '.01', '.03', '.1', '.3'};
        coneModulationsResponseRange = maxConeModulationResponses*[-1 1];
        coneModulationsResponseRangeLog = [0.001 1.05] * maxConeModulationResponses;
    else
        coneModulationsResponseTicks = -1:0.5:1;
        coneModulationsResponseTicksDouble = -1:0.25:1;
        coneModulationsResponseTicksLog = [0.01 0.03 0.1 0.3 1];
        coneModulationsResponseTickLabels = strrep(strrep(sprintf('%+.1f\n', coneModulationsResponseTicks), '0.', '.'), '1.0', '1');
        coneModulationsResponseTickLabelsDouble = strrep(strrep(sprintf('%+.1f\n', coneModulationsResponseTicksDouble), '0.', '.'), '1.0', '1');
        coneModulationsResponseTickLogLabels = {'.01', '.03', '.1', '.3', '1.0'};
        coneModulationsResponseRange = 1*[-1 1];
        coneModulationsResponseRangeLog = [0.001 1.05] * maxConeModulationResponses;
    end
    coneModulationsResponseTickLabels = strrep(coneModulationsResponseTickLabels, '+.0', '0');

    
    if (~mRGCsOperateOnBackgroundAdaptedPhotocurrents)
        photocurrentsReponseRange = maxPhotocurrentResponses;
        photocurrentsResponseTicks = -200:5:200;
        photocurrentsResponseTicksLog = [0.03 0.1 0.3 1 3];
        photocurrentsResponseTickLabels = sprintf('%.0f\n', photocurrentsResponseTicks);
        photocurrentsResponseTickLogLabels = strrep(strrep(sprintf('%.2f\n', photocurrentsResponseTicksLog), '0.', '.'), '00', '0');
        photocurrentsReponseRangeLog = [0.001 1.05]*9;
    else
        if (maxPhotocurrentResponses <= 1)
            photocurrentsResponseTicks = -1:0.5:1;
            photocurrentsResponseTicksDouble = -1:0.25:1;
            photocurrentsResponseTicksLog = [0.1 0.3 1 3];
            photocurrentsResponseTickLabels = strrep(sprintf('%+.1f\n', photocurrentsResponseTicks), '0.', '.');
            photocurrentsResponseTickLabelsDouble = strrep(sprintf('%+.1f\n', photocurrentsResponseTicksDouble), '0.', '.');
            photocurrentsResponseTickLogLabels = sprintf('%.1f\n', photocurrentsResponseTicksLog);
            photocurrentsReponseRange = maxPhotocurrentResponses * [-1 1]*1.05;
            photocurrentsReponseRangeLog = [0.001 1.05]*maxPhotocurrentResponses;
        elseif (maxPhotocurrentResponses <= 2)
            photocurrentsResponseTicks = -2:1:2;
            photocurrentsResponseTicksDouble = -2:0.5:2;
            photocurrentsResponseTicksLog = [0.1 0.3 1 3];
            photocurrentsResponseTickLabels = strrep(sprintf('%+.1f\n', photocurrentsResponseTicks), '0.', '.');
            photocurrentsResponseTickLabelsDouble = strrep(sprintf('%+.1f\n', photocurrentsResponseTicksDouble), '0.', '.');
            photocurrentsResponseTickLogLabels = sprintf('%.1f\n', photocurrentsResponseTicksLog);
            photocurrentsReponseRange = maxPhotocurrentResponses * [-1 1]*1.05;
            photocurrentsReponseRangeLog = [0.001 1.05]*maxPhotocurrentResponses;
        elseif (maxPhotocurrentResponses <= 3)
            photocurrentsResponseTicks = -3:1:3;
            photocurrentsResponseTicksDouble = -3:0.5:3;
            photocurrentsResponseTicksLog = [0.1 0.3 1 3];
            photocurrentsResponseTickLabels = strrep(sprintf('%+.1f\n', photocurrentsResponseTicks), '0.', '.');
            photocurrentsResponseTickLabelsDouble = strrep(sprintf('%+.1f\n', photocurrentsResponseTicksDouble), '0.', '.');
            photocurrentsResponseTickLogLabels = sprintf('%.1f\n', photocurrentsResponseTicksLog);
            photocurrentsReponseRange = maxPhotocurrentResponses * [-1 1]*1.05;
            photocurrentsReponseRangeLog = [0.001 1.05]*maxPhotocurrentResponses;
        elseif (maxPhotocurrentResponses <= 6)
            photocurrentsResponseTicks = -6:2:6;
            photocurrentsResponseTicksDouble = -6:1:6;
            photocurrentsResponseTicksLog = [0.1 0.3 1 3];
            photocurrentsResponseTickLabels = sprintf('%+.0f\n', photocurrentsResponseTicks);
            photocurrentsResponseTickLabelsDouble = sprintf('%+.0f\n', photocurrentsResponseTicksDouble);
            photocurrentsResponseTickLogLabels = sprintf('%.1f\n', photocurrentsResponseTicksLog);
            photocurrentsReponseRange = maxPhotocurrentResponses * [-1 1]*1.05;
            photocurrentsReponseRangeLog = [0.001 1.05]*maxPhotocurrentResponses;
        elseif (maxPhotocurrentResponses <= 8)
            photocurrentsResponseTicks = -8:4:8;
            photocurrentsResponseTicksLog = [0.1 0.3 1 3];
            photocurrentsResponseTickLabels = sprintf('%+.0f\n', photocurrentsResponseTicks);
            photocurrentsResponseTickLogLabels = sprintf('%.1f\n', photocurrentsResponseTicksLog);
            photocurrentsReponseRange = maxPhotocurrentResponses * [-1 1]*1.05;
            photocurrentsReponseRangeLog = [0.001 1.05]*maxPhotocurrentResponses;
        elseif (maxPhotocurrentResponses <= 15)
            photocurrentsResponseTicks = -15:5:15;
            photocurrentsResponseTicksLog = [0.1 0.3 1 3 10];
            photocurrentsResponseTickLabels= sprintf('%+.0f\n', photocurrentsResponseTicks);
            photocurrentsResponseTickLogLabels = sprintf('%.1f\n', photocurrentsResponseTicksLog);
            photocurrentsReponseRange = maxPhotocurrentResponses * [-1 1]*1.05;
            photocurrentsReponseRangeLog = [0.001 1.05]*maxPhotocurrentResponses;
        elseif (maxPhotocurrentResponses <= 30)
            photocurrentsResponseTicks= -30:15:30;
            photocurrentsResponseTicksLog = [0.1 0.3 1 3 10 30];
            photocurrentsResponseTickLabels = sprintf('%+.0f\n', photocurrentsResponseTicks);
            photocurrentsResponseTickLogLabels = sprintf('%.1f\n', photocurrentsResponseTicksLog);
            photocurrentsReponseRange = maxPhotocurrentResponses * [-1 1]*1.05;
            photocurrentsReponseRangeLog = [0.001 1.05]*maxPhotocurrentResponses;
        else
            photocurrentsResponseTicks = -100:20:100;
            photocurrentsResponseTicksLog = [1 3 10 30 100];
            photocurrentsResponseTickLabels = sprintf('%.0f\n', photocurrentsResponseTicks);
            photocurrentsResponseTickLogLabels = sprintf('%.0f\n', photocurrentsResponseTicksLog);
            photocurrentsReponseRange = 100 * [-1 1]*1.05;
            photocurrentsReponseRangeLog = [0.001 1.05]*maxPhotocurrentResponses;
        end
    end
end


function [xx, sfTicks, sfTickLabels] = generateSFticks(sfSupport)
    xx = -numel(sfSupport):1:numel(sfSupport);
    sfTicks = zeros(1, numel(xx));
    sfTickLabels = cell(1, numel(xx));

    for ix = 1:numel(xx)
        sfTicks(ix) = xx(ix);
        if (xx(ix) == 0)
            sfTickLabels{ix} = sprintf('0');
        elseif (xx(ix)>=1)
            theSF = sfSupport(xx(ix));
            if (mod(xx(ix),2) == 0)
                if (theSF >= 5)
                    sfTickLabels{ix} = sprintf('%1.0f', theSF);
                elseif (theSF >= 1)
                    sfTickLabels{ix} = sprintf('%1.1f', theSF);
                elseif (theSF >= 0.1)
                    sfTickLabels{ix} = sprintf('%.1f', theSF);
                elseif (theSF >= 0.01)
                    sfTickLabels{ix} = sprintf('%.2f', theSF);
                else
                    sfTickLabels{ix} = sprintf('%.3f', theSF);
                end
            else
                sfTickLabels{ix} = '';
            end
        else
            theSF = -sfSupport(abs(xx(ix)));
            if (mod(xx(ix),2) == 0)
                if (theSF <= -5)
                    sfTickLabels{ix} = sprintf('%1.0f', theSF);
                elseif (theSF <= -1)
                    sfTickLabels{ix} = sprintf('%1.1f', theSF);
                elseif (theSF <=-0.1)
                    sfTickLabels{ix} = sprintf('%.1f', theSF);
                elseif (theSF <= -0.01)
                    sfTickLabels{ix} = sprintf('%.2f', theSF);
                else
                    sfTickLabels{ix} = sprintf('%.3f', theSF);
                end
            else
                sfTickLabels{ix} = '';
            end
        end

    end

end


function [theAlignedPhotocurrentsResponses, theAlignedConeModulationsResponses, theAlignedTemporalSupportSeconds] = ...
    timeAlignConeModulationAndPhotocurrentResponses(...
       theConeModulationsBasedResponses, thePhotocurrentsBasedResponses, ...
       theConeModulationsBasedResponseTemporalSupportSeconds, ...
       thePhotocurrentsBasedResponseTemporalSupportSeconds, ...
       theConeModulationsBasedSTFphaseSpectra, thePhotocurrentsBasedSTFphaseSpectra, ...
       stimParams, coneResponseExtraSamplesNum, photocurrentResponseExtraSamplesNum, debugPhotocurrentPhaseAlignment)

    % Zero phase align, and temporal alignment of cone modulation and
    % photocurrent responses so we can compare them time point by time point
    theAlignedPhotocurrentsResponses = zeros(numel(stimParams.orientationDegs), numel(stimParams.spatialFrequencyCPD), numel(thePhotocurrentsBasedResponseTemporalSupportSeconds)-photocurrentResponseExtraSamplesNum);
    theAlignedConeModulationsResponses = theAlignedPhotocurrentsResponses;

    for iORI = 1:numel(stimParams.orientationDegs)
        for iSF = 1:numel(stimParams.spatialFrequencyCPD)

            theConeModulationsBasedResponse = squeeze(theConeModulationsBasedResponses(iORI, iSF, :));
            idx1 = 1:(numel(theConeModulationsBasedResponse)-coneResponseExtraSamplesNum);

            thePhotocurrentBasedResponse = squeeze(thePhotocurrentsBasedResponses(iORI, iSF, :));
            idx2 = 1:(numel(thePhotocurrentBasedResponse)-photocurrentResponseExtraSamplesNum);

            % Zero phase alignment for cone modulations
            phaseForAlignment = theConeModulationsBasedSTFphaseSpectra(iORI, iSF);
            thePhaseAlignedConeModulationsBasedResponse = phaseAlignResponse(theConeModulationsBasedResponse(idx1),...
                phaseForAlignment, ...
                theConeModulationsBasedResponseTemporalSupportSeconds(idx1), ...
                1./(stimParams.temporalFrequencyHz), ...
                false);

            % Zero phase alignment for photocurrents
            phaseForAlignment = thePhotocurrentsBasedSTFphaseSpectra(iORI, iSF);
            thePhaseAlignedPhotocurrentBasedResponse = phaseAlignResponse(thePhotocurrentBasedResponse(idx2),...
                phaseForAlignment, ...
                thePhotocurrentsBasedResponseTemporalSupportSeconds(idx2), ...
                1./(stimParams.temporalFrequencyHz), ...
                debugPhotocurrentPhaseAlignment);

            % Interpolate cone modulations to same time base as photocurrent responses
            interpolationMethod = 'linear';
            thePhaseAlignedConeModulationsBasedResponse = interp1(...
                theConeModulationsBasedResponseTemporalSupportSeconds(idx1), ...
                thePhaseAlignedConeModulationsBasedResponse, ...
                thePhotocurrentsBasedResponseTemporalSupportSeconds(idx2), interpolationMethod);

            % Accumulate responses
            theAlignedTemporalSupportSeconds = thePhotocurrentsBasedResponseTemporalSupportSeconds(idx2);
            theAlignedConeModulationsResponses(iORI,iSF,:) = thePhaseAlignedConeModulationsBasedResponse;
            theAlignedPhotocurrentsResponses(iORI,iSF,:) = thePhaseAlignedPhotocurrentBasedResponse;
        end
    end
end


function theResponse = phaseAlignResponse(theUnshiftedResponse, theResponsePhaseDegs, theTemporalSupportSeconds, theStimulusPeriodSeconds, demo)

    dt = theTemporalSupportSeconds(2)-theTemporalSupportSeconds(1);
    fullPeriods = (theTemporalSupportSeconds(end)-theTemporalSupportSeconds(1)+dt)/theStimulusPeriodSeconds;
    
    sizeResponse = size(theUnshiftedResponse);

    theSampleDegs = 360*fullPeriods/(numel(theTemporalSupportSeconds));

    if (demo)
        theResponsePhaseDegs
    end


    if (fullPeriods < 1.5)
        if (theResponsePhaseDegs > 180)
            theResponsePhaseDegs = theResponsePhaseDegs-360;
        end
    else
        %if (theResponsePhaseDegs > 0)
        %    theResponsePhaseDegs = theResponsePhaseDegs-360;
        %end
    end


    if (demo)
        theResponsePhaseDegs
        pause
    end

    
    theShiftAmountSamples = sign(theResponsePhaseDegs) * round(abs(theResponsePhaseDegs)/theSampleDegs);
    theResponse  = circshift(theUnshiftedResponse, -theShiftAmountSamples);

    if (demo)
        figure(1000); clf;
        plot(theTemporalSupportSeconds, theUnshiftedResponse, 'ks');
        hold on;
        for iDemo = 1:abs(theShiftAmountSamples)
            theResponse  = circshift(theUnshiftedResponse, -sign(theShiftAmountSamples)*iDemo);
            plot(theTemporalSupportSeconds, theResponse, 'r-');
            title(sprintf('shift: %d samples', -iDemo))
            pause
        end
        
    end

    theResponse = reshape(theResponse, sizeResponse);
end




