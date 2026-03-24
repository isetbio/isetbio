%
% RGCMosaicAnalyzer.visualize.coneExcitationsVsPhotocurrentsSTF()
%

function coneExcitationsVsPhotocurrentsSTF(...
    theConeModulationsBasedSTFamplitudeSpectra, ...
    thePhotocurrentsBasedSTFamplitudeSpectra, ...
    theConeModulationsBasedSTFphaseSpectra, ...
    thePhotocurrentsBasedSTFphaseSpectra, ...
    theConeModulationsBasedResponses, ...
    thePhotocurrentsBasedResponses, ...
    theConeModulationsBasedResponseTemporalSupportSeconds, ...
    thePhotocurrentsBasedResponseTemporalSupportSeconds, ...
    stimParams, theRGCindex, ...
    maxConeModulationResponses, ...
    maxPhotocurrentResponses, ...
    varargin)


    p = inputParser;
    p.addParameter('exportPDFdirectory', '', @(x)(isempty(x)||ischar(x)));
    % Execute the parser
    p.parse(varargin{:});
    exportPDFdirectory = p.Results.exportPDFdirectory;

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


    hFig = figure(1); clf;
    set(hFig, 'Position', [10 10 2250 850], 'Color', [1 1 1]);


    % Cone modulation based mRGC response time-series
    ax = subplot(2,5,[1 6]);
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
    xlabel('time (seconds)');
    ylabel('mRGC response');
    title(ax, sprintf('cone modulations-based mRGC responses\n(zero phase)'));



    % Photocurrent based mRGC response time-series
    ax = subplot(2,5,[2 7]);
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
    xlabel('time (seconds)');
    ylabel('mRGC response');
    title(ax, sprintf('photocurrents-based mRGC responses\n(phase aligned with cone modulations)'));


    % Cone-based STFs for the examined orientations
    ax = subplot(2,5,4);
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
    grid(ax, 'on');
    legend(ax, theLegends, 'Location', 'SouthWest');
    set(ax, 'FontSize', 16)
    xlabel('spatial frequency (c/deg)');
    ylabel('STF');
    title(ax, 'cone modulations - based STFs');


    % Photocurrents-based STFs for the examined orientations
    ax = subplot(2,5,5);
    hold(ax,'on')
    for iORI = 1:numel(stimParams.orientationDegs)
        plot(stimParams.spatialFrequencyCPD, squeeze(thePhotocurrentsBasedSTFamplitudeSpectra(iORI,:)), ...
            'o-', 'Color', squeeze(oriColors(iORI,:)), 'LineWidth', 1.5, ...
            'MarkerSize', 12, 'MarkerFaceColor', squeeze(oriColors(iORI,:)), 'MarkerEdgeColor', 0.5*squeeze(oriColors(iORI,:)));
    end
    set(ax, 'XScale', 'log', 'XTick', [0.01 0.03 0.1 0.3 1 3 10 30]);
    grid(ax, 'on');
    legend(ax, theLegends, 'Location', 'SouthWest');
    set(ax, 'FontSize', 16)
    xlabel('spatial frequency (c/deg)');
    title(ax, 'photocurrents - based STFs');


    % 2D STF for cone modulations
    ax = subplot(2,5,9);
    imagesc(ax, xx,yy,coneModulationsSTF2D);
    set(ax, 'XTick', [sfTicks(1) 0 sfTicks(end)], 'XTickLabel', {sfTickLabels{1}, '0', sfTickLabels{end}}, ...
        'YTick', sfTicks, 'YTickLabel', sfTickLabels);
    axis(ax, 'image');
    set(ax, 'CLim', [0 max(coneModulationsSTFMatrix(:))]);
    set(ax, 'FontSize', 16)
    title(ax, 'cone modulations-based 2D STF');
    xlabel('spatial frequency (c/deg)');
    ylabel('spatial frequency (c/deg)');

    % 2D STF for photocurrents
    ax = subplot(2,5,10);
    imagesc(ax, xx,yy,photocurrentsSTF2D);
    set(ax, 'XTick', [sfTicks(1) 0 sfTicks(end)], 'XTickLabel', {sfTickLabels{1}, '0', sfTickLabels{end}}, ...
        'YTick', sfTicks', 'YTickLabel', sfTickLabels);
    axis(ax, 'image');
    set(ax, 'CLim', [0 max(photocurrentsSTFMatrix(:))]); 
    set(ax, 'FontSize', 16)
    title(ax, 'photocurrents-based 2D STF')
    xlabel('spatial frequency (c/deg)');
    

    % Plotocurrent - based mRGC response time series (zero phase)
    ax = subplot(2,5,3);
    theYLims = max(abs(thePhotocurrentsBasedResponses(:))) * [-1 1];

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
    xlabel('time (seconds)');
    ylabel('mRGC response');
    title(ax, sprintf('photocurrents-based mRGC responses\n(zero phase)'));


    ax = subplot(2,5,8);
    idx = ~isnan(theConeModulationsBasedResponse);
    allConeModulationsResponses = allConeModulationsResponses(idx);
    allPhotocurrentResponses = allPhotocurrentResponses(idx);
    theXLims = maxConeModulationResponses * [-1 1];
    theYLims = maxPhotocurrentResponses * [-1 1];
    plot(ax, allConeModulationsResponses(:), allPhotocurrentResponses(:), 'k.');
    grid(ax, 'on')
    axis(ax, 'square')
    if (maxPhotocurrentResponses <= 20)
        yTicks = -20:5:20;
    elseif (maxPhotocurrentResponses <= 40)
        yTicks = -40:10:40;
    else
        yTicks = -100:20:100;
    end

    set(ax, 'FontSize', 16, 'Ylim', theYLims, 'XLim', theXLims, 'XTick', -1:0.2:1, 'YTick', yTicks);
    xlabel('cone modulations based');
    ylabel('photocurrents based');
    title(ax, 'photocurrent non-linearity');


    drawnow;

    if (~isempty(exportPDFdirectory))
        theFilename = sprintf('RGC_%d_nominalC_%2.0f%%_%2.0fCDM2_%2.1fHz.pdf', ...
            theRGCindex, stimParams.contrast*100, stimParams.backgroundLuminanceCdM2, stimParams.temporalFrequencyHz);

        % Generate figure dir if it does not exist
        theScriptName = strrep(exportPDFdirectory, ISETBioPaperAndGrantCodeRootDirectory, '');
        theScriptName = strrep(theScriptName, '/','');
        theFiguresDir = ISETBioPaperAndGrantCodeFigureDirForScript(theScriptName);

        NicePlot.exportFigToPDF(fullfile(theFiguresDir,theFilename), hFig, 300);
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
