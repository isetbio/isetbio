function csfPDF(app)

    if (isfield(app.products, 'csfData'))
        hFig = figure('Visible','off'); clf;
        set(hFig, 'Position', [10 10 500 700], 'Color', [1 1 1]);
        ax = subplot('Position', [0.13 0.08 0.83 0.89]);

        % Plot Watson's pyramid of visibility for data of Watson 1987
        [sfSupport, WatsonPyramidOfVisibilitySensitivity] = ...
            CSFGeneratorApp.generate.WatsonPyramidOfVisibilityData(...
            app.products.csfData.spatialFrequencySupport(1), ...
            app.products.csfData.spatialFrequencySupport(end), ...
            app.stimParams.meanLuminanceCdM2, ...
            app.csfParams.constantParameter);
        
        plot(ax, sfSupport, WatsonPyramidOfVisibilitySensitivity, 'k--', 'LineWidth', 1.5);
        hold(ax, 'on');
        plot(ax, app.products.csfData.spatialFrequencySupport, app.products.csfData.sensitivity, 'bo-', ...
                 'MarkerFaceColor', [0.5 0.5 1.0], 'LineWidth', 1.5, ...
                 'MarkerSize', 14); 
             
        lHandle = legend(ax, {...
                        sprintf('Watson''s foveal PoV (%s)', app.csfParams.constantParameter)
                        sprintf('ISETbio (%s @ %2.1f degs, stimSize: %2.1f degs)', app.csfParams.sourceSignal, app.roiParams.radialEccentricityDegs, app.stimParams.sizeDegs)
                        ...
                        }, 'location', 'NorthOutside');
        set(lHandle,'Box','off');
        
        set(ax, 'XLim', [1 40], ...
                'XTick', [1 3 10 30 60 100], 'YTick', [1 3 10 30 100 300 1000 3000 10000], ...
                'XScale', 'log', 'YScale', 'log', 'FontSize', 14);
        set(ax, 'YLim', [1 max([1000 (round(app.products.csfData.sensitivity/1000)+1)*1000])]);
        box(ax, 'on');
        grid(ax, 'on');
        xlabel(ax, 'spatial frequency (c/deg)');
        ylabel(ax, 'contrast sensitivity');

        % Export to PDF
        pdfFileName = sprintf('csf_%d.pdf', app.pdfIndex);
        app.pdfIndex = app.pdfIndex + 1;
        NicePlot.exportFigToPDF(pdfFileName, hFig, 300);
        close(hFig);
        
        % Open generated PDF
        open(pdfFileName);
        
    end
    
    
end
