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
        % Plot the computed CSF
         plot(app.csfView, csfData.spatialFrequencySupport, csfData.sensitivity, 'k-', ...
             'Color', [0.3 0.3 0.3], 'LineWidth', 2);
         hold(app.csfView, 'on');
         for iSF = 1:numel(csfData.spatialFrequencySupport)
             app.csfDataPointHandles(iSF) = scatter(app.csfView, ...
                csfData.spatialFrequencySupport(iSF),  csfData.sensitivity(iSF), 14*14, 'o',  ...
                'MarkerEdgeColor', [0.3 0.3 0.3], ...
                'MarkerFaceColor', squeeze(app.csfLineColors(iSF,:)), ...
                'MarkerFaceAlpha', 0.9, 'LineWidth', 2, ...
                'HandleVisibility','off' ... % do not show legends
                );
         end    
         
         % Plot Watson's pyramid of visibility for data of Watson 1987
         if (strcmp(app.csfParams.constantParameter, 'constant cycles'))
            temporalFrequency = 0.0; cW = 0;
            logLuminanceNits = log10(app.stimParams.meanLuminanceCdM2);
            % Table 1, of Watson 2018, "The Field of View, the Field of Resolution and the Field of Contrast Sensitivity" (Luminance, CCG)
            cF = -0.091;
            cL = 0.391;
            c0 = 1.380;
            sfSupportHiRes = 3:1:60;
            logS = c0 + cW*temporalFrequency + cF*sfSupportHiRes + cL*logLuminanceNits;
            S = 10.^logS;
            plot(app.csfView, sfSupportHiRes, S, 'k--', 'LineWidth', 1.5);
         end
         legend(app.csfView, {'isetbio', 'Watson 2018'});
    end
    
end

function initializeCSFView(app)
    for iSF = 1:numel(app.csfDataPointHandles)
        set(app.csfDataPointHandles(iSF), 'HandleVisibility', 'on');
    end
    cla(app.csfView);
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