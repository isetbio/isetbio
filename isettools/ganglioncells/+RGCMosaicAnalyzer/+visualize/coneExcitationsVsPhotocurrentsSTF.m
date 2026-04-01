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
    p.addParameter('withSuperimposedPSF', [], @(x)(isempty(x)||isstruct(x)));
    p.addParameter('spatialSupportTickSeparationArcMin', 7, @isscalar);
    p.addParameter('extraPlotType', 'differentialPhaseSpectraAndBPIs', @(x)(ismember(x, {'differentialPhaseSpectraAndBPIs', 'individualPhaseSpectra', '2DSTFs'})));
    p.addParameter('exportPDFdirectory', '', @(x)(isempty(x)||ischar(x)));
    p.addParameter('videoOBJ', []);
    % Execute the parser
    p.parse(varargin{:});

    theSuperimposedPSF = p.Results.withSuperimposedPSF;
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
       'leftMargin',     0.04, ...
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
            'withSuperimposedPSF', theSuperimposedPSF, ...
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

    coneModulationsLogResponseTicks = [0.01 0.03 0.1 0.3 1];
    photocurrentsLogResponseTicks = [0.1 0.3 1 3 10 30 100];

    coneModulationsLogResponseTickLabels = strrep(sprintf('%1.2f\n', coneModulationsLogResponseTicks), '0.', '.');
    photocurrentsLogResponseTickLabels = sprintf('%1.1f\n', photocurrentsLogResponseTicks);

    if (maxConeModulationResponses <= 0.2)
        coneModulationsResponseTicks = -0.2:0.05:0.2;
        coneModulationsResponseTicks = coneModulationsResponseTicks(abs(coneModulationsResponseTicks)<=maxConeModulationResponses);
        coneModulationsResponseTickLabels = strrep(sprintf('%+.2f\n', coneModulationsResponseTicks), '0.', '.');
        coneModulationsResponseRange = maxConeModulationResponses*[-1 1];
        
    elseif (maxConeModulationResponses <= 0.4)
        coneModulationsResponseTicks = -0.4:0.1:0.4;
        coneModulationsResponseTicks = coneModulationsResponseTicks(abs(coneModulationsResponseTicks)<=maxConeModulationResponses);
        coneModulationsResponseTickLabels = strrep(sprintf('%+.1f\n', coneModulationsResponseTicks), '0.', '.');
        coneModulationsResponseRange = maxConeModulationResponses*[-1 1];

    elseif (maxConeModulationResponses <= 0.8)
        coneModulationsResponseTicks = -0.8:0.2:0.8;
        coneModulationsResponseTicks = coneModulationsResponseTicks(abs(coneModulationsResponseTicks)<=maxConeModulationResponses);
        coneModulationsResponseTickLabels = strrep(sprintf('%+.1f\n', coneModulationsResponseTicks), '0.', '.');
        coneModulationsResponseRange = maxConeModulationResponses*[-1 1];
    else
        coneModulationsResponseTicks = -1:0.25:1;
        coneModulationsResponseTickLabels = strrep(sprintf('%+.2f\n', coneModulationsResponseTicks), '0.', '.');
        coneModulationsResponseRange = 1*[-1 1];
    end
    coneModulationsResponseTickLabels = strrep(coneModulationsResponseTickLabels, '+.0', '0');

    if (maxPhotocurrentResponses <= 2.5)
        photocurrentsResponseTicks = -2.5:0.5:2.5;
        photocurrentsResponseTickLabels = strrep(sprintf('%.1f\n', photocurrentsResponseTicks), '0.', '.');
        photocurrentsReponseRange = maxPhotocurrentResponses * [-1 1];
    elseif (maxPhotocurrentResponses <= 4)
        photocurrentsResponseTicks = -4:1:4;
        photocurrentsResponseTickLabels = sprintf('%.0f\n', photocurrentsResponseTicks);
        photocurrentsReponseRange = maxPhotocurrentResponses * [-1 1];
    elseif (maxPhotocurrentResponses <= 8)
        photocurrentsResponseTicks = -8:2:8;
        photocurrentsResponseTickLabels = sprintf('%.0f\n', photocurrentsResponseTicks);
        photocurrentsReponseRange = maxPhotocurrentResponses * [-1 1];
    elseif (maxPhotocurrentResponses <= 16)
        photocurrentsResponseTicks = -16:4:16;
        photocurrentsResponseTickLabels= sprintf('%.0f\n', photocurrentsResponseTicks);
        photocurrentsReponseRange = maxPhotocurrentResponses * [-1 1];
    elseif (maxPhotocurrentResponses <= 32)
        photocurrentsResponseTicks= -32:8:32;
        photocurrentsResponseTickLabels = sprintf('%.0f\n', photocurrentsResponseTicks);
        photocurrentsReponseRange = maxPhotocurrentResponses * [-1 1];
    else
        photocurrentsResponseTicks = -100:20:100;
        photocurrentsResponseTickLabels = sprintf('%.0f\n', photocurrentsResponseTicks);
        photocurrentsReponseRange = 100 * [-1 1];
    end


    coneResponseExtraSamplesNum = 1;

    %stimParams
    %theConeModulationsBasedResponseTemporalSupportSeconds
    %halfConeResponseTimeBin = 0.5*(theConeModulationsBasedResponseTemporalSupportSeconds(2)-theConeModulationsBasedResponseTemporalSupportSeconds(1))
    %[~,photocurrentLastIndexForPhaseAlignment] = min(abs(thePhotocurrentsBasedResponseTemporalSupportSeconds-(theConeModulationsBasedResponseTemporalSupportSeconds(end-coneResponseExtraSamplesNum)+halfConeResponseTimeBin)));
    %photocurrentResponseExtraSamplesNum = numel(thePhotocurrentsBasedResponseTemporalSupportSeconds)-photocurrentLastIndexForPhaseAlignment
    photocurrentResponseExtraSamplesNum = 0;

    % Cone modulation based mRGC response time-series
    ax = subplot('Position', subplotPosVectors(1,2).v);

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
            idx = 1:(numel(theConeModulationsBasedResponse)-coneResponseExtraSamplesNum);
            thePhaseAlignedConeModulationsBasedResponse = phaseAlignResponse(theConeModulationsBasedResponse(idx),...
                phaseForAlignment, ...
                theConeModulationsBasedResponseTemporalSupportSeconds(idx), ...
                1./(stimParams.temporalFrequencyHz), ...
                false);
            stairs(ax, theConeModulationsBasedResponseTemporalSupportSeconds(idx), thePhaseAlignedConeModulationsBasedResponse, ...
                'Color', sfColors(5+iSF,:), 'LineWidth', 1.5);
        end
    end
    grid(ax, 'on')
    axis(ax, 'square')
    set(ax, 'FontSize', 20, 'Ylim', coneModulationsResponseRange, 'YTick', coneModulationsResponseTicks, 'YTickLabel', coneModulationsResponseTickLabels);
    set(ax, 'XLim', [0 theConeModulationsBasedResponseTemporalSupportSeconds(end)]);
    set(ax, 'XTick', timeTicks);
    
    xtickangle(ax, 0)
    xlabel(ax,'time (seconds)',  'FontAngle', 'italic');
    ylabel(ax,'mRGC response (modulation)',  'FontAngle', 'italic');
    title(ax, sprintf('cone modulations-based mRGC responses\n(zero phase)'), 'FontWeight', 'normal', 'FontSize', 16);



   
    % Photocurrent based mRGC response time-series
    ax = subplot('Position', subplotPosVectors(2,2).v);
    
    hold (ax, 'on');
    for iORI = 1:numel(stimParams.orientationDegs)
        if (iORI == 1)
            sfColors = brewermap(5+numel(stimParams.spatialFrequencyCPD), 'reds');
        else
            sfColors = brewermap(5+numel(stimParams.spatialFrequencyCPD), 'blues');
        end
        for iSF = 1:numel(stimParams.spatialFrequencyCPD)
            % Align with respect to the cone modulation response to
            % reveal temporal delay of the photocurrent response
            % with respect to the coneModulation response
            phaseForAlignment = theConeModulationsBasedSTFphaseSpectra(iORI, iSF);
            thePhotocurrentBasedResponse = squeeze(thePhotocurrentsBasedResponses(iORI, iSF, :));
            idx = 1:(numel(thePhotocurrentBasedResponse)-photocurrentResponseExtraSamplesNum);
            thePhaseAlignedPhotocurrentBasedResponse = phaseAlignResponse(thePhotocurrentBasedResponse(idx),...
                phaseForAlignment, ...
                thePhotocurrentsBasedResponseTemporalSupportSeconds(idx), ...
                1./(stimParams.temporalFrequencyHz), ...
                false);
            plot(ax, thePhotocurrentsBasedResponseTemporalSupportSeconds(idx), thePhaseAlignedPhotocurrentBasedResponse, ...
                '-', 'Color', sfColors(5+iSF,:), 'LineWidth', 1.5);
        end
    end
    grid(ax, 'on')
    axis(ax, 'square')
    set(ax, 'FontSize', 20, 'Ylim', photocurrentsReponseRange,  'XLim', [0 theConeModulationsBasedResponseTemporalSupportSeconds(end)]);
    set(ax, 'XTick', timeTicks, 'YTick', photocurrentsResponseTicks);
    xtickangle(ax, 0)
    xlabel(ax,'time (seconds)',  'FontAngle', 'italic');
    ylabel(ax,'mRGC response (pAmps)',  'FontAngle', 'italic');
    title(ax, sprintf('photocurrents-based mRGC responses\n(phase relative to cone modulations response)'), 'FontWeight', 'normal', 'FontSize', 16);


    % Plotocurrent - based mRGC response time series (zero phase)
    ax = subplot('Position', subplotPosVectors(1,3).v);

    allPhotocurrentsResponses = zeros(numel(stimParams.orientationDegs), numel(stimParams.spatialFrequencyCPD), numel(thePhotocurrentsBasedResponseTemporalSupportSeconds)-photocurrentResponseExtraSamplesNum);
    allConeModulationsResponses = allPhotocurrentsResponses;

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
            idx1 = 1:(numel(theConeModulationsBasedResponse)-coneResponseExtraSamplesNum);
            thePhaseAlignedConeModulationsBasedResponse = phaseAlignResponse(theConeModulationsBasedResponse(idx1),...
                phaseForAlignment, ...
                theConeModulationsBasedResponseTemporalSupportSeconds(idx1), ...
                1./(stimParams.temporalFrequencyHz), ...
                false);

            phaseForAlignment = thePhotocurrentsBasedSTFphaseSpectra(iORI, iSF);
            thePhotocurrentBasedResponse = squeeze(thePhotocurrentsBasedResponses(iORI, iSF, :));
            idx2 = 1:(numel(thePhotocurrentBasedResponse)-photocurrentResponseExtraSamplesNum);
            thePhaseAlignedPhotocurrentBasedResponse = phaseAlignResponse(thePhotocurrentBasedResponse(idx2),...
                phaseForAlignment, ...
                thePhotocurrentsBasedResponseTemporalSupportSeconds(idx2), ...
                1./(stimParams.temporalFrequencyHz), ...
                false);
            plot(ax, thePhotocurrentsBasedResponseTemporalSupportSeconds(idx2), thePhaseAlignedPhotocurrentBasedResponse, ...
                '-', 'Color', sfColors(5+iSF,:), 'LineWidth', 1.5);

            % Interpolate cone modulations to same time base as photocurrent responses
            interpolationMethod = 'linear';
            thePhaseAlignedConeModulationsBasedResponse = interp1(...
                theConeModulationsBasedResponseTemporalSupportSeconds(idx1), ...
                thePhaseAlignedConeModulationsBasedResponse, ...
                thePhotocurrentsBasedResponseTemporalSupportSeconds(idx2), interpolationMethod);

            % Accumulate responses
            allConeModulationsResponses(iORI,iSF,:) = thePhaseAlignedConeModulationsBasedResponse;
            allPhotocurrentsResponses(iORI,iSF,:) = thePhaseAlignedPhotocurrentBasedResponse;
            
        end
    end
    grid(ax, 'on')
    axis(ax, 'square')
    set(ax, 'FontSize', 20, 'Ylim', photocurrentsReponseRange, 'XLim', [0 theConeModulationsBasedResponseTemporalSupportSeconds(end)]);
    xtickangle(ax, 0)
    set(ax, 'XTick', timeTicks, 'YTick', photocurrentsResponseTicks);
    xlabel(ax,'time (seconds)',  'FontAngle', 'italic');
    ylabel(ax,'mRGC response (pAmps)',  'FontAngle', 'italic');
    title(ax, sprintf('photocurrents-based mRGC responses\n(zero phase)'), 'FontWeight', 'normal', 'FontSize', 16);


    % The NONLINEAR relationship between cone excitations and photocurrents
    ax = subplot('Position', subplotPosVectors(2,3).v);

    lineAlpha = 0.8;
    lineWidth = 1.5;

    for iORI = 1:numel(stimParams.orientationDegs)
        if (iORI == 1)
            lineColors = brewermap(5+numel(stimParams.spatialFrequencyCPD), 'reds');
        else
            lineColors = brewermap(5+numel(stimParams.spatialFrequencyCPD), 'blues');
        end
        for iSF = 1:numel(stimParams.spatialFrequencyCPD)
            theConeModulationsResponseTimeSeries = allConeModulationsResponses(iORI,iSF,:);
            thePhotocurrentsResponseTimeSeries = allPhotocurrentsResponses(iORI,iSF,:);
    
            p = patch(ax, [theConeModulationsResponseTimeSeries(:);NaN],[thePhotocurrentsResponseTimeSeries(:);NaN],'k');
            set(p, 'edgeColor', lineColors(5+iSF,:), 'edgeAlpha', lineAlpha, 'lineWidth', lineWidth);
        end
    end

    hold(ax, 'on');

    plot(ax, [0 0], photocurrentsReponseRange, 'k-', 'LineWidth', 1.0);
    plot(ax, coneModulationsResponseRange, [0 0], 'k-', 'LineWidth', 1.0);
    hold(ax, 'off')
    grid(ax, 'on')
    axis(ax, 'square');
    set(ax, 'FontSize', 20, ...
        'Ylim', photocurrentsReponseRange, 'XLim', coneModulationsResponseRange, ...
        'XTick', coneModulationsResponseTicks, 'XTickLabel', coneModulationsResponseTickLabels, ...
        'YTick', photocurrentsResponseTicks, 'YTickLabel', photocurrentsResponseTickLabels);
    xtickangle(ax, 0)
    xlabel(ax,sprintf('cone modulations-based\nmRGC responses (modulation)'),  'FontAngle', 'italic');
    ylabel(ax,sprintf('photocurrents-based\nmRGC responses (pAmps)'),  'FontAngle', 'italic');
    title(ax, sprintf('photocurrent non-linearity\n(%2.2fHz, %2.0fcd/m2, %2.0f%%)', stimParams.temporalFrequencyHz, stimParams.backgroundLuminanceCdM2, 100*stimParams.contrast'), ...
        'FontWeight', 'normal', 'FontSize', 16);



    % Cone-based STF amplitude spectra for the examined orientations
    ax = subplot('Position', subplotPosVectors(1,4).v);
    hold(ax,'on')

    thePlotHandles = [];
    for iORI = 1:numel(stimParams.orientationDegs)
        if (iORI == 1)
            lineColors = brewermap(5+numel(stimParams.spatialFrequencyCPD), 'reds');
        else
            lineColors = brewermap(5+numel(stimParams.spatialFrequencyCPD), 'blues');
        end
        plot(ax,stimParams.spatialFrequencyCPD, squeeze(theConeModulationsBasedSTFamplitudeSpectra(iORI,:)), ...
                '-', 'Color', squeeze(lineColors(5+7,:)), 'LineWidth', 1.5);

        for iSF = 1:numel(stimParams.spatialFrequencyCPD)
            p = plot(ax,stimParams.spatialFrequencyCPD(iSF), squeeze(theConeModulationsBasedSTFamplitudeSpectra(iORI,iSF)), ...
                'o-', 'Color', squeeze(lineColors(5+iSF,:)), 'LineWidth', 1.5, ...
                'MarkerSize', 12, 'MarkerFaceColor', squeeze(lineColors(5+iSF,:)), 'MarkerEdgeColor', 0.5*squeeze(lineColors(5+iSF,:)));
            if (iSF == 7)
                thePlotHandles(iORI) = p;
                theLegends{iORI} = sprintf('%d degs (BPI:%2.2f)', stimParams.orientationDegs(iORI), thePhotocurrentsBasedBPIs(end, iORI));
            end
        end
        
    end
    axis(ax, 'square')
    set(ax, 'XScale', 'log', 'XTick', [0.01 0.03 0.1 0.3 1 3 10 30]);
    set(ax, 'YScale', 'log', 'YTick', coneModulationsLogResponseTicks, 'YTickLabel', coneModulationsLogResponseTickLabels);
    set(ax, 'XLim', [0.05 60], 'YLim', [0.01 coneModulationsResponseRange(2)]);
    grid(ax, 'on');
    legend(ax, thePlotHandles, theLegends, 'Location', 'SouthWest');
    set(ax, 'FontSize', 20)
    xtickangle(ax, 0)
    xlabel(ax,'spatial frequency (c/deg)',  'FontAngle', 'italic');
    ylabel(ax,'STF amplitude (modulation)',  'FontAngle', 'italic');
    title(ax, sprintf('cone modulations - based\nSTF amplitude spectra'), 'FontWeight', 'normal', 'FontSize', 16);
    

    % Photocurrents-based STF amplitude spectra for the examined orientations
    ax = subplot('Position', subplotPosVectors(1,5).v);
    hold(ax,'on')

    thePlotHandles = [];
    for iORI = 1:numel(stimParams.orientationDegs)
        if (iORI == 1)
            lineColors = brewermap(5+numel(stimParams.spatialFrequencyCPD), 'reds');
        else
            lineColors = brewermap(5+numel(stimParams.spatialFrequencyCPD), 'blues');
        end
        plot(ax,stimParams.spatialFrequencyCPD, squeeze(thePhotocurrentsBasedSTFamplitudeSpectra(iORI,:)), ...
                '-', 'Color', squeeze(lineColors(5+7,:)), 'LineWidth', 1.5);

        for iSF = 1:numel(stimParams.spatialFrequencyCPD)
            p = plot(ax,stimParams.spatialFrequencyCPD(iSF), squeeze(thePhotocurrentsBasedSTFamplitudeSpectra(iORI,iSF)), ...
                'o-', 'Color', squeeze(lineColors(5+iSF,:)), 'LineWidth', 1.5, ...
                'MarkerSize', 12, 'MarkerFaceColor', squeeze(lineColors(5+iSF,:)), 'MarkerEdgeColor', 0.5*squeeze(lineColors(5+iSF,:)));
            if (iSF == 7)
                thePlotHandles(iORI) = p;
                theLegends{iORI} = sprintf('%d degs (BPI:%2.2f)', stimParams.orientationDegs(iORI), thePhotocurrentsBasedBPIs(end, iORI));
            end
        end
        
    end
    axis(ax, 'square')
    
    set(ax, 'YScale', 'log', 'YTick', photocurrentsLogResponseTicks, 'YTickLabel', photocurrentsLogResponseTickLabels);
    set(ax, 'XScale', 'log', 'XTick', [0.01 0.03 0.1 0.3 1 3 10 30]);
    set(ax, 'XLim', [0.05 60], 'YLim', [0.1 photocurrentsReponseRange(2)]);
    grid(ax, 'on');
    legend(ax, thePlotHandles, theLegends, 'Location', 'SouthWest');
    set(ax, 'FontSize', 20)
    xtickangle(ax, 0)
    xlabel(ax,'spatial frequency (c/deg)',  'FontAngle', 'italic');
    ylabel(ax,'STF amplitude (pAmps)',  'FontAngle', 'italic');
    title(ax, sprintf('photocurrents - based\nSTF amplitude spectra'), 'FontWeight', 'normal', 'FontSize', 16);


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
            set(ax, 'XLim', [0.05 60], 'YLim', 120+60*[-1 1]);
            grid(ax, 'on');
            legend(ax, thePlotHandles, theLegends, 'Location', 'SouthEAST');
            set(ax, 'FontSize', 20)
            xtickangle(ax, 0)
            xlabel(ax,'spatial frequency (c/deg)',  'FontAngle', 'italic');
            ylabel(ax,'STF phase difference (degs)',  'FontAngle', 'italic');
            title(ax,sprintf('differential STF phase spectra\n(photocurrents - cone modulations)'), 'FontWeight', 'normal', 'FontSize', 16);


            % The BPI's for all cells up to this one
            ax = subplot('Position', subplotPosVectors(2,5).v);
            hold(ax,'on')
            theLegends = cell(1,numel(stimParams.orientationDegs));
        
            oriColors = brewermap(numel(stimParams.orientationDegs), 'RdBu');

            pHandles = [];

            % All previous cell BPIs (small disks)
            for iORI = 1:numel(stimParams.orientationDegs)
                scatter(ax,theConeModulationsBasedBPIs(1:end-1, iORI), thePhotocurrentsBasedBPIs(1:end-1, iORI), 11^2, ...
                    'MarkerFaceColor', squeeze(oriColors(iORI,:)), 'MarkerFaceAlpha', 0.5, ...
                    'MarkerEdgeColor', squeeze(oriColors(iORI,:)), 'LineWidth', 1.5, 'MarkerEdgeAlpha', 0.9);
            end

            % Black line joining the BPIs for all orientations for the current cell
            plot(ax,theConeModulationsBasedBPIs(end, :), thePhotocurrentsBasedBPIs(end, :), 'k-', 'LineWidth', 3.0);
            plot(ax,theConeModulationsBasedBPIs(end, :), thePhotocurrentsBasedBPIs(end, :), '-','Color', [0.1 0.8 1.0], 'LineWidth', 1.0);

            % Current cell BPIs (large cyan square)
            for iORI = 1:numel(stimParams.orientationDegs)
                pHandles(iORI) = scatter(ax,theConeModulationsBasedBPIs(end, iORI), thePhotocurrentsBasedBPIs(end, iORI), 20^2, ...
                    'Marker', 's', 'MarkerFaceColor', [0.1 0.8 1.0] , 'MarkerFaceAlpha', 1.0, ...
                    'MarkerEdgeColor', squeeze(oriColors(iORI,:)), 'MarkerEdgeAlpha', 1.0, 'LineWidth', 2.0);

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
            title(ax, sprintf('STF bandpass index\n(STF(0)/STF(peak sf)'), 'FontWeight', 'normal', 'FontSize', 16);

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
            title(ax, 'cone modulations-based 2D STF', 'FontWeight', 'normal', 'FontSize', 16);
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
            title(ax, 'photocurrents-based 2D STF',  'FontWeight', 'normal', 'FontSize', 16)
          

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

function theResponse = phaseAlignResponse(theUnshiftedResponse, theResponsePhaseDegs, theTemporalSupportSeconds, theStimulusPeriodSeconds, demo)

    dt = theTemporalSupportSeconds(2)-theTemporalSupportSeconds(1);
    fullPeriods = (theTemporalSupportSeconds(end)-theTemporalSupportSeconds(1)+dt)/theStimulusPeriodSeconds;
    
    sizeResponse = size(theUnshiftedResponse);

    theSampleDegs = 360*fullPeriods/(numel(theTemporalSupportSeconds));

    if (demo)
        theResponsePhaseDegs
    end


    if (theStimulusPeriodSeconds < 1.5)
        if (theResponsePhaseDegs > 180)
            theResponsePhaseDegs = theResponsePhaseDegs-360;
        end
    else
        if (theResponsePhaseDegs > 0)
            theResponsePhaseDegs = theResponsePhaseDegs-360;
        end
    end


    if (demo)
        theResponsePhaseDegs
        pause
    end

    

    theShiftAmountSamples = sign(theResponsePhaseDegs) * round(abs(theResponsePhaseDegs)/theSampleDegs);
    theResponse  = circshift(theUnshiftedResponse, -theShiftAmountSamples);

    if (demo)
        theShiftAmountSamples
        figure(1000); clf;
        plot(theTemporalSupportSeconds, theUnshiftedResponse, 'ks');
        hold on;
        for iDemo = 1:abs(theShiftAmountSamples)
            theResponse  = circshift(theUnshiftedResponse, -iDemo);
            plot(theTemporalSupportSeconds, theResponse, 'r-');
            title(sprintf('shift: %d samples', -iDemo))
            pause
        end
        
    end

    theResponse = reshape(theResponse, sizeResponse);
end




