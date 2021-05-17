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
        for iSF = 1:numel(sfSupport)
            scatter(app.csfView, sfSupport(iSF), 1.2, 14*14, 'v', ...
                     'LineWidth', 1.5, 'MarkerFaceColor', squeeze(app.csfLineColors(iSF,:)), ...
                     'MarkerFaceAlpha', 0.3, 'MarkerEdgeColor', [0.2 0.2 0.2], ...
                     'MarkerEdgeAlpha', 0.3);
            hold(app.csfView, 'on');
        end
    else
         plot(app.csfView, csfData.spatialFrequencySupport, csfData.sensitivity, 'k-', ...
             'Color', [0.3 0.3 0.3], 'LineWidth', 2);
         hold(app.csfView, 'on');
         for iSF = 1:numel(csfData.spatialFrequencySupport)
             scatter(app.csfView, ...
                csfData.spatialFrequencySupport(iSF),  csfData.sensitivity(iSF), 14*14, 'o',  ...
                'MarkerEdgeColor', [0.3 0.3 0.3], ...
                'MarkerFaceColor', squeeze(app.csfLineColors(iSF,:)), ...
                'MarkerFaceAlpha', 0.9, 'LineWidth', 2);
         end    
    end
    
end

function initializeCSFView(app)
    set(app.csfView, 'XLim', [0.5 100], 'YLim', [1 1000], ...
                'XTick', [1 3 10 30 60 100], 'YTick', [1 3 10 30 100 300 1000 3000 10000], ...
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