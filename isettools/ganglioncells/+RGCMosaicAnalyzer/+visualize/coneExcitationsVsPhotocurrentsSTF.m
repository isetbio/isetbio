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
    varargin)


    p = inputParser;
    p.addParameter('exportPDFdirectory', '', @(x)(isempty(x)||ischar(x)));
    p.addParameter('videoOBJ', []);
    % Execute the parser
    p.parse(varargin{:});
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

    spatialSupportTickSeparationArcMin = 6;

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

   
    % Add a scale bar, 0.093 degs in size
    % which is around 25 microns at an eccentricity of 25 degs
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
    set(ax, 'FontSize', 16);
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
            'yLimsRange', [-0.7 1.01], ...
            'axesToRenderIn', ax);
    set(ax, 'FontSize', 16);
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
    set(ax, 'FontSize', 16, 'Ylim', theYLims, 'XLim', [0 theConeModulationsBasedResponseTemporalSupportSeconds(end)]);
    xtickangle(ax, 0)
    xlabel(ax,'time (seconds)');
    ylabel(ax,'mRGC response');
    title(ax, sprintf('cone modulations-based mRGC responses\n(zero phase)'));


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
    set(ax, 'FontSize', 16, 'Ylim', theYLims, 'XLim', [0 theConeModulationsBasedResponseTemporalSupportSeconds(end)]);
    xtickangle(ax, 0)
    xlabel(ax,'time (seconds)');
    ylabel(ax,'mRGC response');
    title(ax, sprintf('photocurrents-based mRGC responses\n(relative phase to cone modulations response)'));


    % Plotocurrent - based mRGC response time series (zero phase)
    ax = subplot('Position', subplotPosVectors(1,3).v);
    theYLims = maxPhotocurrentResponses * [-1 1];

    allConeModulationResponses = [];
    allPhotocurrentResponses = [];

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

            % Find time correspondence between cone modulations and
            % photocurrents
            interpolationMethod = 'nearest';
            theConeModulationsBasedResponse = interp1(theConeModulationsBasedResponseTemporalSupportSeconds, theConeModulationsBasedResponse, ...
                thePhotocurrentsBasedResponseTemporalSupportSeconds, interpolationMethod);
            allPhotocurrentResponses(numel(allPhotocurrentResponses)+(1:numel(thePhotocurrentBasedResponse))) = thePhotocurrentBasedResponse;
            allConeModulationsResponses(numel(allConeModulationResponses)+(1:numel(theConeModulationsBasedResponse))) = theConeModulationsBasedResponse;
        end
    end
    grid(ax, 'on')
    set(ax, 'FontSize', 16, 'Ylim', theYLims, 'XLim', [0 theConeModulationsBasedResponseTemporalSupportSeconds(end)]);
    xtickangle(ax, 0)
    xlabel(ax,'time (seconds)');
    ylabel(ax,'mRGC response');
    title(ax, sprintf('photocurrents-based mRGC responses\n(zero phase)'));


    % The relationship between cone excitations and photocurrents
    ax = subplot('Position', subplotPosVectors(2,3).v);
    idx = ~isnan(theConeModulationsBasedResponse);
    allConeModulationsResponses = allConeModulationsResponses(idx);
    allPhotocurrentResponses = allPhotocurrentResponses(idx);
    
    if (maxPhotocurrentResponses <= 2.5)
        yTicks = -2.5:0.5:2.5;
        theYLims = 2.5 * [-1 1];
    elseif (maxPhotocurrentResponses <= 5)
        yTicks = -5:1:5;
        theYLims = 5 * [-1 1];
    elseif (maxPhotocurrentResponses <= 10)
        yTicks = -10:2:10;
        theYLims = 10 * [-1 1];
    elseif (maxPhotocurrentResponses <= 20)
        yTicks = -20:5:20;
        theYLims = 20 * [-1 1];
    elseif (maxPhotocurrentResponses <= 40)
        yTicks = -40:10:40;
        theYLims = 10 * [-1 1];
    else
        yTicks = -100:20:100;
        theYLims = 100 * [-1 1];
    end

    if (maxConeModulationResponses <= 0.2)
        xTicks = -0.2:0.05:0.2;
        theXLims = 0.2*[-1 1];
    elseif (maxConeModulationResponses <= 0.5)
        xTicks = -0.5:0.1:0.5;
        theXLims = 0.5*[-1 1];
    elseif (maxConeModulationResponses <= 0.8)
        xTicks = -0.8:0.2:0.8;
        theXLims = 0.8*[-1 1];
    else
        xTicks = -1:0.25:1;
        theXLims = 1*[-1 1];
    end

    scatter(ax, allConeModulationsResponses(:), allPhotocurrentResponses(:), 81, ...
        'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 0.5, 'MarkerEdgeColor', [1 0 0], 'MarkerFaceColor', [1 0.5 0.5]);
    hold(ax, 'on');
    plot(ax, [0 0], theYLims, 'k-', 'LineWidth', 1.0);
    plot(ax, theXLims, [0 0], 'k-', 'LineWidth', 1.0);
    hold(ax, 'off')
    grid(ax, 'on')
    axis(ax, 'square');
    set(ax, 'FontSize', 16, 'Ylim', theYLims, 'XLim', theXLims, 'XTick', xTicks, 'YTick', yTicks);
    xtickangle(ax, 0)
    xlabel(ax,'cone modulations-based mRGC responses');
    ylabel(ax,'photocurrents-based mRGC responses');
    title(ax, 'photocurrent non-linearity');



    % Cone-based STFs for the examined orientations
    ax = subplot('Position', subplotPosVectors(1,4).v);
    hold(ax,'on')
    theLegends = cell(1,numel(stimParams.orientationDegs));

    oriColors = brewermap(numel(stimParams.orientationDegs), 'RdBu');

    for iORI = 1:numel(stimParams.orientationDegs)
        plot(stimParams.spatialFrequencyCPD, squeeze(theConeModulationsBasedSTFamplitudeSpectra(iORI,:)), ...
            'o-', 'Color', squeeze(oriColors(iORI,:)), 'LineWidth', 1.5, ...
            'MarkerSize', 12, 'MarkerFaceColor', squeeze(oriColors(iORI,:)), 'MarkerEdgeColor', 0.5*squeeze(oriColors(iORI,:)));
        theLegends{iORI} = sprintf('%d degs', stimParams.orientationDegs(iORI));
    end
    set(ax, 'XScale', 'log', 'XTick', [0.01 0.03 0.1 0.3 1 3 10 30]);
    set(ax, 'XLim', [0.01 60], 'YLim', [0 maxConeModulationResponses]);
    grid(ax, 'on');
    legend(ax, theLegends, 'Location', 'SouthWest');
    set(ax, 'FontSize', 16)
    xtickangle(ax, 0)
    xlabel(ax,'spatial frequency (c/deg)');
    ylabel(ax,'STF');
    title(ax, 'cone modulations - based STFs');


    % Photocurrents-based STFs for the examined orientations
    ax = subplot('Position', subplotPosVectors(1,5).v);
    hold(ax,'on')
    for iORI = 1:numel(stimParams.orientationDegs)
        plot(stimParams.spatialFrequencyCPD, squeeze(thePhotocurrentsBasedSTFamplitudeSpectra(iORI,:)), ...
            'o-', 'Color', squeeze(oriColors(iORI,:)), 'LineWidth', 1.5, ...
            'MarkerSize', 12, 'MarkerFaceColor', squeeze(oriColors(iORI,:)), 'MarkerEdgeColor', 0.5*squeeze(oriColors(iORI,:)));
    end
    set(ax, 'XScale', 'log', 'XTick', [0.01 0.03 0.1 0.3 1 3 10 30]);
    set(ax, 'XLim', [0.01 60], 'YLim', [0 maxPhotocurrentResponses]);
    grid(ax, 'on');
    legend(ax, theLegends, 'Location', 'SouthWest');
    set(ax, 'FontSize', 16)
    xtickangle(ax, 0)
    xlabel(ax,'spatial frequency (c/deg)');
    title(ax, 'photocurrents - based STFs');


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
    set(ax, 'FontSize', 16)
    xtickangle(ax, 0)
    title(ax, 'cone modulations-based 2D STF');
    xlabel(ax,'spatial frequency (c/deg)');
    ylabel(ax,'spatial frequency (c/deg)');

    

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
    set(ax, 'FontSize', 16)
    xtickangle(ax, 0)
    title(ax, 'photocurrents-based 2D STF')
    xlabel(ax,'spatial frequency (c/deg)');
    

    
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
