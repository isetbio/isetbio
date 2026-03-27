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
    stimParams, theRGCindex, theMRGCmosaic, ...
    maxConeModulationResponses, ...
    maxPhotocurrentResponses, ...
    theConeModulationsBasedBPIs, ...
    thePhotocurrentsBasedBPIs, ...
    varargin)


    p = inputParser;
    p.addParameter('spatialSupportTickSeparationArcMin', 7, @isscalar);
    p.addParameter('extraPlotType', 'differentialPhaseSpectraAndBPIs', @(x)(ismember(x, {'differentialPhaseSpectraAndBPIs', 'individualPhaseSpectra', '2DSTFs'})));
    p.addParameter('exportPDFdirectory', '', @(x)(isempty(x)||ischar(x)));
    p.addParameter('videoOBJ', []);
    % Execute the parser
    p.parse(varargin{:});


    spatialSupportTickSeparationArcMin = p.Results.spatialSupportTickSeparationArcMin;
    extraPlotType = p.Results.extraPlotType;
    exportPDFdirectory = p.Results.exportPDFdirectory;
    videoOBJ = p.Results.videoOBJ;

    x = [];
    y = [];
    coneModulationsSTFMatrix = [];
    photocurrentsSTFMatrix = [];


    for iSF = 1:numel(stimParams.spatialFrequencyCPD)
        for iORI = 1:numel(stimParams.orientationDegs)
            x(numel(x)+1) = iSF * cosd(stimParams.orientationDegs(iORI));
            y(numel(y)+1) = iSF * sind(stimParams.orientationDegs(iORI));
            coneModulationsSTFMatrix(numel(coneModulationsSTFMatrix)+1) = theConeModulationsBasedSTFamplitudeSpectra(iORI,iSF);
            photocurrentsSTFMatrix(numel(photocurrentsSTFMatrix)+1) = thePhotocurrentsBasedSTFamplitudeSpectra(iORI, iSF);
            x(numel(x)+1) = -x(numel(x));
            y(numel(y)+1) = -y(numel(y));
            coneModulationsSTFMatrix(numel(coneModulationsSTFMatrix)+1) = theConeModulationsBasedSTFamplitudeSpectra(iORI,iSF);
            photocurrentsSTFMatrix(numel(photocurrentsSTFMatrix)+1) = thePhotocurrentsBasedSTFamplitudeSpectra(iORI, iSF);
        end
    end

    interpolationMethod = 'natural';
    extrapolationMethod = 'none';
    FconeModulations = scatteredInterpolant(x(:), y(:), coneModulationsSTFMatrix(:), interpolationMethod, extrapolationMethod);
    Fphotocurrents = scatteredInterpolant(x(:), y(:), photocurrentsSTFMatrix(:), interpolationMethod, extrapolationMethod);

    [xx, sfTicks, sfTickLabels] = generateSFticks(stimParams.spatialFrequencyCPD);

    yy = xx;
    [X,Y] = meshgrid(xx,yy);

    coneModulationsSTF2D = FconeModulations(X(:),Y(:));
    coneModulationsSTF2D = reshape(coneModulationsSTF2D, numel(xx), numel(yy));

    photocurrentsSTF2D = Fphotocurrents(X(:),Y(:));
    photocurrentsSTF2D = reshape(photocurrentsSTF2D, numel(xx), numel(yy));

    subplotPosVectors = NicePlot.getSubPlotPosVectors(...
       'rowsNum', 2, ...
       'colsNum', 5, ...
       'heightMargin',  0.11, ...
       'widthMargin',    0.06, ...
       'leftMargin',     0.03, ...
       'rightMargin',    0.00, ...
       'bottomMargin',   0.07, ...
       'topMargin',      0.05);


    hFig = figure(1); clf;
    set(hFig, 'Position', [10 10 2240 960], 'Color', [1 1 1]);

    if (~isempty(exportPDFdirectory))
        thePDFFilename = sprintf('RGC_%d_nominalC_%2.0f%%_%2.0fCDM2_%2.1fHz.pdf', ...
            theRGCindex, stimParams.contrast*100, stimParams.backgroundLuminanceCdM2, stimParams.temporalFrequencyHz);

        % Generate figure dir if it does not exist
        theScriptName = strrep(exportPDFdirectory, ISETBioPaperAndGrantCodeRootDirectory, '');
        theScriptName = strrep(theScriptName, '/','');
        theFiguresDir = ISETBioPaperAndGrantCodeFigureDirForScript(theScriptName);
    end

    % The cone pooling map
    ax = subplot('Position', subplotPosVectors(1,1).v);

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
            'axesToRenderIn', ax, ...
            'minConeWeightForVisualizingRFcenterPooling', minCenterConeWeightForConePoolingMapVisualization, ...
            'minConeWeightForVisualizingRFsurroundPooling', minSurroundConeWeightForConePoolingMapVisualization, ...
            'minSurroundConeWeightRelativity', 'center', ...
            'withLineWeightingFunctions', false, ...
            'domainVisualizationLimits', domainVisualizationLimitsSingleRFDefault, ...
            'domainVisualizationTicks', domainVisualizationTicksSingleRFDefault, ...
            'spatialSupportTickSeparationArcMin', spatialSupportTickSeparationArcMin, ...
            'scaleBarDegs', scaleBarDegs, ...
            'doNotLabelScaleBar', true, ...
            'plotTitle', sprintf('RGC #%d', theRGCindex), ...
            'figNo', [], ...
            'figPos', [], ...
            'exportVisualizationPDF',false, ...
            'exportVisualizationPNG',false);
    set(ax, 'FontSize', 20);
    xtickangle(ax, 0)


    % The horizontal profile of the cone pooling map
    ax = subplot('Position', subplotPosVectors(2,1).v);
    whichMeridian = 'horizontal';		% choose between {'horizontal', 'vertical'}
    spatialSupportCenterDegs = theMRGCmosaic.rgcRFpositionsDegs(theRGCindex,:);
	RGCMosaicConstructor.visualize.centerAndSurroundConePoolingLineWeightingFunctions(...
            '', [], [], ...
		    spatialSupportCenterDegs, spatialSupportTickSeparationArcMin, ...
		    centerLineWeightingFunctions, surroundLineWeightingFunctions, whichMeridian, ...
            'domainVisualizationLimits', domainVisualizationLimitsSingleRFDefault(1:2), ...
            'domainVisualizationTicks', domainVisualizationTicksSingleRFDefault.x, ...
            'yLimsRange', [-0.8 1.01], ...
            'axesToRenderIn', ax);
    set(ax, 'FontSize', 20);
    xtickangle(ax, 0)

    % Cone modulation based mRGC response time-series
    ax = subplot('Position', subplotPosVectors(1,2).v);
    theYLims = maxConeModulationResponses * [-1 1];
    hold (ax, 'on');
    for iORI = 1:numel(stimParams.orientationDegs)
        if (iORI == 1)
            sfColors = brewermap(5+numel(stimParams.spatialFrequencyCPD), 'reds');
        else
            sfColors = brewermap(5+numel(stimParams.spatialFrequencyCPD), 'blues');
        end
        for iSF = 1:numel(stimParams.spatialFrequencyCPD)
            theConeModulationsBasedResponse = squeeze(theConeModulationsBasedResponses(iORI, iSF, :));
            phaseForAlignment = theConeModulationsBasedSTFphaseSpectra(iORI, iSF);
            theConeModulationsBasedResponse = phaseAlignResponse(theConeModulationsBasedResponse,...
                phaseForAlignment, ...
                theConeModulationsBasedResponseTemporalSupportSeconds, ...
                1);
            stairs(ax, theConeModulationsBasedResponseTemporalSupportSeconds, theConeModulationsBasedResponse, ...
                'Color', sfColors(5+iSF,:), 'LineWidth', 1.5);
        end
    end
    grid(ax, 'on')
    axis(ax, 'square')
    set(ax, 'FontSize', 20, 'Ylim', theYLims, 'XLim', [0 theConeModulationsBasedResponseTemporalSupportSeconds(end)]);
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
    set(ax, 'XTick', timeTicks);
    yTicks = -1:0.1:1;
    set(ax, 'YTick', yTicks, 'YTickLabel', sprintf('%.1f\n', yTicks));
    xtickangle(ax, 0)
    xlabel(ax,'time (seconds)',  'FontAngle', 'italic');
    ylabel(ax,'mRGC response',  'FontAngle', 'italic');
    title(ax, sprintf('cone modulations-based mRGC responses\n(zero phase)'), 'FontWeight', 'normal');


    % Photocurrent based mRGC response time-series
    ax = subplot('Position', subplotPosVectors(2,2).v);
    theYLims = maxPhotocurrentResponses * [-1 1];
    hold (ax, 'on');
    for iORI = 1:numel(stimParams.orientationDegs)
        if (iORI == 1)
            sfColors = brewermap(5+numel(stimParams.spatialFrequencyCPD), 'reds');
        else
            sfColors = brewermap(5+numel(stimParams.spatialFrequencyCPD), 'blues');
        end
        for iSF = 1:numel(stimParams.spatialFrequencyCPD)
            phaseForAlignment = thePhotocurrentsBasedSTFphaseSpectra(iORI, iSF);

            % Align with respect to the cone modulation response to
            % reveal temporal delay of the photocurrent response
            % with respect to the coneModulation response
            phaseForAlignment = theConeModulationsBasedSTFphaseSpectra(iORI, iSF);

            thePhotocurrentBasedResponse = squeeze(thePhotocurrentsBasedResponses(iORI, iSF, :));
            thePhotocurrentBasedResponse = phaseAlignResponse(thePhotocurrentBasedResponse,...
                phaseForAlignment, ...
                thePhotocurrentsBasedResponseTemporalSupportSeconds, ...
                1);
            plot(ax, thePhotocurrentsBasedResponseTemporalSupportSeconds, thePhotocurrentBasedResponse, ...
                '-', 'Color', sfColors(5+iSF,:), 'LineWidth', 1.5);
        end
    end
    grid(ax, 'on')
    axis(ax, 'square')
    set(ax, 'FontSize', 20, 'Ylim', theYLims, 'XLim', [0 theConeModulationsBasedResponseTemporalSupportSeconds(end)]);
    set(ax, 'XTick', timeTicks);
    xtickangle(ax, 0)
    xlabel(ax,'time (seconds)',  'FontAngle', 'italic');
    ylabel(ax,'mRGC response',  'FontAngle', 'italic');
    title(ax, sprintf('photocurrents-based mRGC responses\n(phase relative to cone modulations response)'), 'FontWeight', 'normal');


    % Plotocurrent - based mRGC response time series (zero phase)
    ax = subplot('Position', subplotPosVectors(1,3).v);
    theYLims = maxPhotocurrentResponses * [-1 1];

    allConeModulationsResponses = [];
    allPhotocurrentsResponses = [];

    hold (ax, 'on');
    for iORI = 1:numel(stimParams.orientationDegs)
        if (iORI == 1)
            sfColors = brewermap(5+numel(stimParams.spatialFrequencyCPD), 'reds');
        else
            sfColors = brewermap(5+numel(stimParams.spatialFrequencyCPD), 'blues');
        end
        for iSF = 1:numel(stimParams.spatialFrequencyCPD)
            phaseForAlignment = theConeModulationsBasedSTFphaseSpectra(iORI, iSF);
            theConeModulationsBasedResponse = squeeze(theConeModulationsBasedResponses(iORI, iSF, :));
            theConeModulationsBasedResponse = phaseAlignResponse(theConeModulationsBasedResponse,...
                phaseForAlignment, ...
                theConeModulationsBasedResponseTemporalSupportSeconds, ...
                1);

            phaseForAlignment = thePhotocurrentsBasedSTFphaseSpectra(iORI, iSF);
            thePhotocurrentBasedResponse = squeeze(thePhotocurrentsBasedResponses(iORI, iSF, :));
            thePhotocurrentBasedResponse = phaseAlignResponse(thePhotocurrentBasedResponse,...
                phaseForAlignment, ...
                thePhotocurrentsBasedResponseTemporalSupportSeconds, ...
                0);
            plot(ax, thePhotocurrentsBasedResponseTemporalSupportSeconds, thePhotocurrentBasedResponse, ...
                '-', 'Color', sfColors(5+iSF,:), 'LineWidth', 1.5);

            % Interpolate cone modulations to same time base as photocurrent responses
            interpolationMethod = 'linear';
            theConeModulationsBasedResponse = interp1(theConeModulationsBasedResponseTemporalSupportSeconds, theConeModulationsBasedResponse, ...
                thePhotocurrentsBasedResponseTemporalSupportSeconds, interpolationMethod);

            % Accumulate responses
            allPhotocurrentsResponses(size(allPhotocurrentsResponses,1)+1,:) = thePhotocurrentBasedResponse;
            allConeModulationsResponses(size(allConeModulationsResponses,1)+1,:) = theConeModulationsBasedResponse;
        end
    end
    grid(ax, 'on')
    axis(ax, 'square')
    set(ax, 'FontSize', 20, 'Ylim', theYLims, 'XLim', [0 theConeModulationsBasedResponseTemporalSupportSeconds(end)]);
    xtickangle(ax, 0)
    set(ax, 'XTick', timeTicks);
    xlabel(ax,'time (seconds)',  'FontAngle', 'italic');
    ylabel(ax,'mRGC response',  'FontAngle', 'italic');
    title(ax, sprintf('photocurrents-based mRGC responses\n(zero phase)'), 'FontWeight', 'normal');


    % The relationship between cone excitations and photocurrents
    ax = subplot('Position', subplotPosVectors(2,3).v);
    %idx = ~isnan(allConeModulationsResponses);
    %allConeModulationsResponses = allConeModulationsResponses(idx);
    %allPhotocurrentsResponses = allPhotocurrentsResponses(idx);
    
    if (maxPhotocurrentResponses <= 2.5)
        yTicks = -2.5:0.5:2.5;
        yTickLabels = strrep(sprintf('%.1f\n', yTicks), '0.', '.');
        theYLims = maxPhotocurrentResponses * [-1 1];
    elseif (maxPhotocurrentResponses <= 4)
        yTicks = -4:1:4;
        yTickLabels = sprintf('%.0f\n', yTicks);
        theYLims = maxPhotocurrentResponses * [-1 1];
    elseif (maxPhotocurrentResponses <= 8)
        yTicks = -8:2:8;
        yTickLabels = sprintf('%.0f\n', yTicks);
        theYLims = maxPhotocurrentResponses * [-1 1];
    elseif (maxPhotocurrentResponses <= 16)
        yTicks = -16:4:16;
        yTickLabels = sprintf('%.0f\n', yTicks);
        theYLims = maxPhotocurrentResponses * [-1 1];
    elseif (maxPhotocurrentResponses <= 32)
        yTicks = -32:8:32;
        yTickLabels = sprintf('%.0f\n', yTicks);
        theYLims = maxPhotocurrentResponses * [-1 1];
    else
        yTicks = -100:20:100;
        yTickLabels = sprintf('%.0f\n', yTicks);
        theYLims = 100 * [-1 1];
    end

    if (maxConeModulationResponses <= 0.2)
        xTicks = -0.2:0.05:0.2;
        xTicks = xTicks(abs(xTicks)<=maxConeModulationResponses);
        xTickLabels = strrep(sprintf('%+.2f\n', xTicks), '0.', '.');
        theXLims = maxConeModulationResponses*[-1 1];
        
    elseif (maxConeModulationResponses <= 0.4)
        xTicks = -0.4:0.1:0.4;
        xTicks = xTicks(abs(xTicks)<=maxConeModulationResponses);
        xTickLabels = strrep(sprintf('%+.1f\n', xTicks), '0.', '.');
        theXLims = maxConeModulationResponses*[-1 1];
        
    elseif (maxConeModulationResponses <= 0.6)
        xTicks = -0.6:0.15:0.6;
        xTicks = xTicks(abs(xTicks)<=maxConeModulationResponses);
        xTickLabels = strrep(sprintf('%+.2f\n', xTicks), '0.', '.');
        theXLims = maxConeModulationResponses*[-1 1];

    elseif (maxConeModulationResponses <= 0.8)
        xTicks = -0.8:0.2:0.8;
        xTicks = xTicks(abs(xTicks)<=maxConeModulationResponses);
        xTickLabels = strrep(sprintf('%+.1f\n', xTicks), '0.', '.');
        theXLims = maxConeModulationResponses*[-1 1];
    else
        xTicks = -1:0.25:1;
        xTickLabels = strrep(sprintf('%+.2f\n', xTicks), '0.', '.');
        theXLims = 1*[-1 1];
    end

    xTickLabels = strrep(xTickLabels, '+.0', '0');

    %scatter(ax, allConeModulationsResponses(:), allPhotocurrentsResponses(:), 64, ...
    %    'MarkerFaceAlpha', 0.1, 'MarkerEdgeAlpha', 0.0, 'MarkerEdgeColor', [1 0 0], 'MarkerFaceColor', [1 0.5 0.5]);

    lineColor = [1 0.35 0.0];
    lineAlpha = 0.2;
    lineWidth = 1.5;
    for iStim = 1:size(allPhotocurrentsResponses,1)
        xCoords = allConeModulationsResponses(iStim,:);
        yCoords = allPhotocurrentsResponses(iStim,:);

        p = patch(ax, [xCoords(:);NaN],[yCoords(:);NaN],'k');
        set(p, 'edgeColor', lineColor, 'edgeAlpha', lineAlpha, 'lineWidth', lineWidth);
    end

    hold(ax, 'on');

    plot(ax, [0 0], theYLims, 'k-', 'LineWidth', 1.0);
    plot(ax, theXLims, [0 0], 'k-', 'LineWidth', 1.0);
    hold(ax, 'off')
    grid(ax, 'on')
    axis(ax, 'square');
    set(ax, 'FontSize', 20, 'Ylim', theYLims, 'XLim', theXLims, ...
        'XTick', xTicks, 'XTickLabel', xTickLabels, ...
        'YTick', yTicks, 'YTickLabel', yTickLabels);
    xtickangle(ax, 0)
    xlabel(ax,sprintf('cone modulations-based\nmRGC responses)'),  'FontAngle', 'italic');
    ylabel(ax,sprintf('photocurrents-based\nmRGC responses'),  'FontAngle', 'italic');
    title(ax, sprintf('photocurrent non-linearity\n(%2.2fHz, %2.0fcd/m2, %2.0f%%)', stimParams.temporalFrequencyHz, stimParams.backgroundLuminanceCdM2, 100*stimParams.contrast'), ...
        'FontWeight', 'normal', 'Color', lineColor);


    % Cone-based STF amplitude spectra for the examined orientations
    ax = subplot('Position', subplotPosVectors(1,4).v);
    hold(ax,'on')
    theLegends = cell(1,numel(stimParams.orientationDegs));

    oriColors = brewermap(numel(stimParams.orientationDegs), 'RdBu');

    for iORI = 1:numel(stimParams.orientationDegs)
        plot(ax, stimParams.spatialFrequencyCPD, squeeze(theConeModulationsBasedSTFamplitudeSpectra(iORI,:)), ...
            'o-', 'Color', squeeze(oriColors(iORI,:)), 'LineWidth', 1.5, ...
            'MarkerSize', 12, 'MarkerFaceColor', squeeze(oriColors(iORI,:)), 'MarkerEdgeColor', 0.5*squeeze(oriColors(iORI,:)));
        theLegends{iORI} = sprintf('%d degs (BPI:%2.2f)', stimParams.orientationDegs(iORI), theConeModulationsBasedBPIs(end, iORI));
    end
    yTicks = 0:0.05:1;
    axis(ax, 'square')
    set(ax, 'XScale', 'log', 'XTick', [0.01 0.03 0.1 0.3 1 3 10 30], 'YTick', yTicks, 'YTickLabel', sprintf('%.2f\n', yTicks));
    set(ax, 'XLim', [0.05 60], 'YLim', [0 maxConeModulationResponses]);
    grid(ax, 'on');
    legend(ax, theLegends, 'Location', 'NorthWest');
    set(ax, 'FontSize', 20)
    xtickangle(ax, 0)
    xlabel(ax,'spatial frequency (c/deg)',  'FontAngle', 'italic');
    ylabel(ax,'STF amplitude',  'FontAngle', 'italic');
    title(ax, sprintf('cone modulations - based\nSTF amplitude spectra'), 'FontWeight', 'normal');


    % Photocurrents-based STF amplitude spectra for the examined orientations
    ax = subplot('Position', subplotPosVectors(1,5).v);
    hold(ax,'on')
    for iORI = 1:numel(stimParams.orientationDegs)
        plot(ax,stimParams.spatialFrequencyCPD, squeeze(thePhotocurrentsBasedSTFamplitudeSpectra(iORI,:)), ...
            'o-', 'Color', squeeze(oriColors(iORI,:)), 'LineWidth', 1.5, ...
            'MarkerSize', 12, 'MarkerFaceColor', squeeze(oriColors(iORI,:)), 'MarkerEdgeColor', 0.5*squeeze(oriColors(iORI,:)));
        theLegends{iORI} = sprintf('%d degs (BPI:%2.2f)', stimParams.orientationDegs(iORI), thePhotocurrentsBasedBPIs(end, iORI));
    end
    axis(ax, 'square')
    set(ax, 'XScale', 'log', 'XTick', [0.01 0.03 0.1 0.3 1 3 10 30]);
    set(ax, 'XLim', [0.05 60], 'YLim', [0 maxPhotocurrentResponses]);
    grid(ax, 'on');
    legend(ax, theLegends, 'Location', 'NorthWest');
    set(ax, 'FontSize', 20)
    xtickangle(ax, 0)
    xlabel(ax,'spatial frequency (c/deg)',  'FontAngle', 'italic');
    title(ax, sprintf('photocurrents - based\nSTF amplitude spectra'), 'FontWeight', 'normal');


    switch (extraPlotType)
        case 'differentialPhaseSpectraAndBPIs'
            ax = subplot('Position', subplotPosVectors(2,4).v);
            hold(ax,'on')
            theLegends = cell(1,numel(stimParams.orientationDegs));
        
            oriColors = brewermap(numel(stimParams.orientationDegs), 'RdBu');

            % Tolerance for a phase jump between consecutive angles that when we exeed,
            % we shift the angles by adding multiples of 360 degs
            toleranceDegs = 180;
            toleranceRadians = toleranceDegs/180*pi;

            for iORI = 1:numel(stimParams.orientationDegs)
                theConeModulationBasedSTFphaseDegs = squeeze(theConeModulationsBasedSTFphaseSpectra(iORI,:));
                thePhotocurrentBasedSTFphaseDegs = squeeze(thePhotocurrentsBasedSTFphaseSpectra(iORI,:));

                % To radians
                theConeModulationBasedSTFphaseRadians = theConeModulationBasedSTFphaseDegs/180*pi;
                thePhotocurrentBasedSTFphaseRadians = thePhotocurrentBasedSTFphaseDegs/180*pi;

                % Dividing complex numbers in Euler form (r * exp(i * theta) involves dividing their magnitudes (r)
                % and subtracting their angles 
                theEulerRatio = exp(1j*thePhotocurrentBasedSTFphaseRadians) ./ exp(1j*theConeModulationBasedSTFphaseRadians);
                thePhaseDifferenceRadians = angle(theEulerRatio);

                % Unwrap and back to degrees
                theUnwrappedPhaseDifferenceDegs = unwrap(thePhaseDifferenceRadians, toleranceRadians)/pi*180;

                plot(ax,stimParams.spatialFrequencyCPD, theUnwrappedPhaseDifferenceDegs, ...
                    'o-', 'Color', squeeze(oriColors(iORI,:)), 'LineWidth', 1.5, ...
                    'MarkerSize', 12, 'MarkerFaceColor', squeeze(oriColors(iORI,:)), 'MarkerEdgeColor', 0.5*squeeze(oriColors(iORI,:)));
                theLegends{iORI} = sprintf('%d degs', stimParams.orientationDegs(iORI));
            end

            axis(ax, 'square')
            yTicks = -360:30:360;
            set(ax, 'XScale', 'log', 'XTick', [0.01 0.03 0.1 0.3 1 3 10 30]);
            set(ax, 'YTick', yTicks);
            set(ax, 'XLim', [0.05 60], 'YLim', [0 270]);
            grid(ax, 'on');
            legend(ax, theLegends, 'Location', 'SouthEast');
            set(ax, 'FontSize', 20)
            xtickangle(ax, 0)
            xlabel(ax,'spatial frequency (c/deg)',  'FontAngle', 'italic');
            ylabel(ax,'STF phase difference (degs)',  'FontAngle', 'italic');
            title(ax,sprintf('differential STF phase spectra\n(photocurrents - cone modulations)'), 'FontWeight', 'normal');


            % The BPI's for all cells up to this one
            ax = subplot('Position', subplotPosVectors(2,5).v);
            hold(ax,'on')
            theLegends = cell(1,numel(stimParams.orientationDegs));
        
            oriColors = brewermap(numel(stimParams.orientationDegs), 'RdBu');

            pHandles = [];
            for iORI = 1:numel(stimParams.orientationDegs)
                scatter(ax,theConeModulationsBasedBPIs(1:end-1, iORI), thePhotocurrentsBasedBPIs(1:end-1, iORI), 11^2, ...
                    'MarkerFaceColor', 0.5*[1 1 1] + 0.5*squeeze(oriColors(iORI,:)), 'MarkerFaceAlpha', 0.7, ...
                    'MarkerEdgeColor', squeeze(oriColors(iORI,:)), 'LineWidth', 1.5, 'MarkerEdgeAlpha', 0.7);
                pHandles(iORI) = scatter(ax,theConeModulationsBasedBPIs(end, iORI), thePhotocurrentsBasedBPIs(end, iORI), 16^2, ...
                    'Marker', 's', 'MarkerFaceColor', 0.5*[1 1 1] + 0.5*squeeze(oriColors(iORI,:)), 'MarkerFaceAlpha', 1.0, ...
                    'MarkerEdgeColor', squeeze(oriColors(iORI,:)), 'MarkerEdgeAlpha', 1.0, 'LineWidth', 1.5);
                theLegends{iORI} = sprintf('%d degs', stimParams.orientationDegs(iORI));
            end
            plot([0 1], [0 1], 'k-', 'LineWidth', 1.5);
            axis(ax, 'square');
            set(ax, 'XTick', 0:0.2:1, 'XTickLabel', {'0.0', '0.2', '0.4', '0.6', '0.8', 'LP'});
            set(ax, 'YTick', 0:0.2:1, 'YTickLabel', {'0.0', '0.2', '0.4', '0.6', '0.8', 'LP'});
            set(ax, 'XLim', [0 1], 'YLim', [0 1]);
            grid(ax, 'on');
            legend(ax, pHandles, theLegends, 'Location', 'SouthEast');
            set(ax, 'FontSize', 20)
            xtickangle(ax, 0)
            xlabel(ax,'cone modulations-based STF', 'FontAngle', 'italic');
            ylabel(ax,'photocurrents-based STF',  'FontAngle', 'italic');
            title(ax, sprintf('STF bandpass indices\n(%d mRGCs)',size(thePhotocurrentsBasedBPIs,1)), 'FontWeight', 'normal');

        case 'individualPhaseSpectra'
            % Cone-based STF phase spectra for the examined orientations
            ax = subplot('Position', subplotPosVectors(2,4).v);
            hold(ax,'on')
            theLegends = cell(1,numel(stimParams.orientationDegs));
        
            oriColors = brewermap(numel(stimParams.orientationDegs), 'RdBu');
        
            % Tolerance for a phase jump between consecutive angles that when we exeed,
            % we shift the angles by adding multiples of 360 degs
            toleranceDegs = 180;
            toleranceRadians = toleranceDegs/180*pi;
    
            for iORI = 1:numel(stimParams.orientationDegs)
                theWrappedPhaseDegs = squeeze(theConeModulationsBasedSTFphaseSpectra(iORI,:));
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
            set(ax, 'FontSize', 20)
            xtickangle(ax, 0)
            xlabel(ax,'spatial frequency (c/deg)',  'FontAngle', 'italic');
            ylabel(ax,'STF phase (degs)',  'FontAngle', 'italic');
        

            % Photocurrents-based STF phase spectra for the examined orientations
            ax = subplot('Position', subplotPosVectors(2,5).v);
            hold(ax,'on')
            for iORI = 1:numel(stimParams.orientationDegs)
                theWrappedPhaseDegs = squeeze(thePhotocurrentsBasedSTFphaseSpectra(iORI,:));
                theUnwrappedPhaseDegs = unwrap(theWrappedPhaseDegs/180*pi, toleranceRadians)/pi*180;
                plot(ax, stimParams.spatialFrequencyCPD, theUnwrappedPhaseDegs, ...
                    'o-', 'Color', squeeze(oriColors(iORI,:)), 'LineWidth', 1.5, ...
                    'MarkerSize', 12, 'MarkerFaceColor', squeeze(oriColors(iORI,:)), 'MarkerEdgeColor', 0.5*squeeze(oriColors(iORI,:)));
            end
            axis(ax, 'square')
            set(ax, 'XScale', 'log', 'XTick', [0.01 0.03 0.1 0.3 1 3 10 30]);
            set(ax, 'XLim', [0.05 60], 'YLim', [-360 360]);
            set(ax, 'YTick', yTicks);
            grid(ax, 'on');
            legend(ax, theLegends, 'Location', 'SouthWest');
            set(ax, 'FontSize', 20)
            xtickangle(ax, 0)
            xlabel(ax,'spatial frequency (c/deg)',  'FontAngle', 'italic');
       

        case '2DSTFs'

            % 2D STF for cone modulations
            ax = subplot('Position', subplotPosVectors(2,4).v);
            %imagesc(ax, xx,yy,coneModulationsSTF2D);
        
            zLevels = (0:0.1:1.0)*max(coneModulationsSTFMatrix(:));
            contourf(ax, xx, yy, coneModulationsSTF2D, zLevels);
        
            set(ax, 'XTick', [sfTicks(1) 0 sfTicks(end)], 'XTickLabel', {sfTickLabels{1}, '0', sfTickLabels{end}}, ...
                'YTick', sfTicks, 'YTickLabel', sfTickLabels);
            colormap(ax, brewermap(1024, '*greys'));
            axis(ax, 'image');
            set(ax, 'CLim', [0 max(coneModulationsSTFMatrix(:))], 'Color', [0 0 0]);
            set(ax, 'FontSize', 20)
            xtickangle(ax, 0)
            title(ax, 'cone modulations-based 2D STF', 'FontWeight', 'normal');
            xlabel(ax,'spatial frequency (c/deg)',  'FontAngle', 'italic');
            ylabel(ax,'spatial frequency (c/deg)',  'FontAngle', 'italic');

    
            % 2D STF for photocurrents
            ax = subplot('Position', subplotPosVectors(2,5).v);
            %imagesc(ax, xx,yy,photocurrentsSTF2D);
            zLevels = (0:0.1:1.0)*max(photocurrentsSTFMatrix(:));
            contourf(ax, xx, yy, photocurrentsSTF2D, zLevels);
        
            set(ax, 'XTick', [sfTicks(1) 0 sfTicks(end)], 'XTickLabel', {sfTickLabels{1}, '0', sfTickLabels{end}}, ...
                'YTick', sfTicks', 'YTickLabel', sfTickLabels);
            colormap(ax, brewermap(1024, '*greys'));
            axis(ax, 'image');
            set(ax, 'CLim', [0 max(photocurrentsSTFMatrix(:))], 'Color', [0 0 0]); 
            set(ax, 'FontSize', 20)
            xtickangle(ax, 0)
            xlabel(ax,'spatial frequency (c/deg)',  'FontAngle', 'italic');
            title(ax, 'photocurrents-based 2D STF',  'FontWeight', 'normal')
          

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

function theResponse = phaseAlignResponse(theResponse, theResponsePhaseDegs, theTemporalSupportSeconds, extraSamples)

    sizeResponse = size(theResponse);
    theSampleDegs = 360/(numel(theTemporalSupportSeconds)-extraSamples);
    if (theResponsePhaseDegs > 180)
        theResponsePhaseDegs = -(360-theResponsePhaseDegs);
    end

    theShiftAmountSamples = sign(theResponsePhaseDegs) * round(abs(theResponsePhaseDegs)/theSampleDegs);

    if (extraSamples > 0)
        theResponse = theResponse(1:end-extraSamples);
    end
    theResponse  = circshift(theResponse , -theShiftAmountSamples);

    if (extraSamples > 0)
        theResponse = [theResponse(:);  theResponse(1:extraSamples)];
    end

    theResponse = reshape(theResponse, sizeResponse);
end




