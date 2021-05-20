function csfView(app, mode, varargin)
    p = inputParser;
    p.addParameter('withData',  []);
    p.parse(varargin{:});
    csfData = p.Results.withData;
            
    switch (mode)
        case 'initialize'
            initializeCSFView(app);
        case 'update'
            updateCSFViewWithNewData(app, csfData);
    end
end

function  updateCSFViewWithNewData(app, csfData)
    hold(app.csfView, 'off');
    % Plot triangles to indicate the spatial frequencies to be examined
    sfSupport = logspace( log10(app.csfParams.spatialFrequencyMin), ...
                          log10(app.csfParams.spatialFrequencyMax), ...
                          app.csfParams.spatialFrequencySamples);
                      
    if (isempty(csfData))
        % Just plot the triangles showing the SF support
        for iSF = 1:numel(sfSupport)
            scatter(app.csfView, sfSupport(iSF), 1.2, 14*14, 'v', ...
                     'LineWidth', 1.5, 'MarkerFaceColor', squeeze(app.csfLineColors(iSF,:)), ...
                     'MarkerFaceAlpha', 0.3, 'MarkerEdgeColor', [0.2 0.2 0.2], ...
                     'MarkerEdgeAlpha', 0.3);
            hold(app.csfView, 'on');
        end
    else
         % Plot Watson's pyramid of visibility for data of Watson 1987
         [sfSupport, WatsonPyramidOfVisibilitySensitivity] = ...
             CSFGeneratorApp.generate.WatsonPyramidOfVisibilityData(...
                    csfData.spatialFrequencySupport(1), ...
                    csfData.spatialFrequencySupport(end), ...
                    app.stimParams.meanLuminanceCdM2, ...
                    app.csfParams.constantParameter);
        
         % Plot Watson's pyramid of visibility
         plot(app.csfView, sfSupport, WatsonPyramidOfVisibilitySensitivity, 'k--', 'LineWidth', 1.5);
         hold(app.csfView, 'on');
         
         % Plot the computed CSF
         plot(app.csfView, csfData.spatialFrequencySupport, csfData.sensitivity, 'k-', ...
             'Color', [0.3 0.3 0.3], 'LineWidth', 2);
         
         
         % Plot the data points
         for iSF = 1:numel(csfData.spatialFrequencySupport)
             app.csfDataPointHandles(iSF) = scatter(app.csfView, ...
                csfData.spatialFrequencySupport(iSF),  csfData.sensitivity(iSF), 14*14, 'o',  ...
                'MarkerEdgeColor', [0.3 0.3 0.3], ...
                'MarkerFaceColor', squeeze(app.csfLineColors(iSF,:)), ...
                'MarkerFaceAlpha', 0.9, 'LineWidth', 2 ...
                );
         end    
         
         set(app.csfView, 'YLim', [1 max([1000 (round(csfData.sensitivity/1000)+1)*1000])]);
         lHandle = legend(app.csfView, {...
                                sprintf('Watson''s foveal PoV (%s)', app.csfParams.constantParameter)...
                                sprintf('ISETbio (%s @ %2.1f degs)', app.csfParams.sourceSignal, app.roiParams.radialEccentricityDegs) ...
             },'Location', 'SouthWest');
         set(lHandle,'Box','off')
    end
    
end

function initializeCSFView(app)
    cla(app.csfView);
    set(app.csfView, 'XLim', [0.5 90], 'YLim', [1 2000], ...
                'XTick', [1 3 10 30 60 100], ...
                'YTick', [1 3 10 30 100 300 1000 3000 10000], ...
                'YTickLabel', {'1', '3', '10', '30', '100', '300', '1k', '3k', '10k'}, ...
                'XScale', 'log', 'YScale', 'log', 'FontSize', 14);
    box(app.csfView, 'on');
    grid(app.csfView, 'on');
    xlabel(app.csfView, 'spatial frequency (c/deg)');
    ylabel(app.csfView, 'contrast sensitivity');

    % Do not show the interactions toolbax
    app.csfView.Toolbar.Visible = 'off';

    % Only allow paning
    app.csfView.Interactions = [];
end