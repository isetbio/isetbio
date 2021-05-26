function opticsView(app, mode)

    switch (mode)
        case 'initialize'
            initializeOpticsView(app);
        case 'update'
            updateOpticsViewWithNewData(app);
    end
end

function initializeOpticsView(app)

    cla(app.opticsView);
    cMap = brewermap(512, '*spectral');
    colormap(app.opticsView, cMap);
    
    app.psfDensityPlotHandle = imagesc(app.opticsView, [-0.01 0.01], [-0.01 0.01], [1 1; 1 1], [0 1]);
    hold(app.opticsView, 'on');
    for k = 1:app.centralConeOutlinesNum
        app.coneOutlineOnPSFPlotHandles(k) = patch(app.opticsView, [0 0],[0 0], [0.5 0.5 0.1], 'FaceAlpha', 0.2, 'EdgeColor', [.3 .3 0.1], 'EdgeAlpha', 1.0, 'LineWidth', 1.0);
    end

    % Plot the cross hairs
    plot(app.opticsView, [0 0], [-20 20], 'k-', 'LineWidth', 1.5, 'Color', [0.5 0.5 0.5]);
    plot(app.opticsView, [-20 20], [0 0], 'k-', 'LineWidth', 1.5, 'Color', [0.5 0.5 0.5]);
    hold(app.opticsView, 'off');   
    box(app.opticsView, 'off');
    title(app.opticsView, 'PSF and cone apertures', 'FontSize', 14, 'FontWeight', 'Normal', 'Color', [0.4 0.4 0.4]);
    
    xlabel(app.opticsView, 'arc min.');
    set(app.opticsView, 'CLim', [0 0.9], 'XColor', [0.3 0.3 0.3], 'YColor', 'none');
    
    set(app.opticsView, ...
        'XLim', 7*[-1 1], ...
        'YLim', 7*[-1 1], ...
        'XTick', -20:2:20, 'YTick', []);
    
    axis(app.opticsView, 'image');
     
    % Do not show the interactions toolbax
    app.opticsView.Toolbar.Visible = 'off';
            
    % Only allow zooming
    app.opticsView.Interactions = [panInteraction zoomInteraction];
            
    % Add listener to zoom-events
    addlistener(app.opticsView, {'XLim', 'YLim'}, 'PostSet', @app.handleOpticsZoomEvent);
    
    
    % opticsViewPSFGrid initialization
    cla(app.opticsViewPSFGrid);
    hold(app.opticsViewPSFGrid, 'on');
    for xo = -1:1
    for yo = -1:1
        % Plot the cross-hairs
        plot(app.opticsViewPSFGrid, xo + 0.4*[-1 1], [yo yo], 'k-', 'LineWidth', 1.0, 'Color', [0.5 0.5 0.5]);
        plot(app.opticsViewPSFGrid, [xo xo], yo + 0.4*[-1 1], 'k-', 'LineWidth', 1.0, 'Color', [0.5 0.5 0.5]);
    end
    end
    hold(app.opticsViewPSFGrid, 'off');
    set(app.opticsViewPSFGrid, 'XLim', [-1.5 1.5], 'YLim', [-1.5 1.5]);
    set(app.opticsViewPSFGrid, 'XColor', 'none', 'YColor', 'none', 'Color', 'none');
    title(app.opticsViewPSFGrid, 'PSF spatial grid', 'FontSize', 14, 'FontWeight', 'Normal', 'Color', [0.4 0.4 0.4]);
    
    
    % Do not show the interactions toolbax
    app.opticsViewPSFGrid.Toolbar.Visible = 'off';
            
    % No interactions
    app.opticsViewPSFGrid.Interactions = [];
end


function updateOpticsViewWithNewData(app)
    % Find index of visualized wavelength
    [~, wIdx] = min(abs(app.components.psf.supportWavelength-app.opticsParams.visualizedWavelength));
    
    % Extract PSF at the desired wavelength
    wavePSF = real(squeeze(app.components.psf.data(:,:,wIdx)));
       
    % Normalize to unit amplitude
    wavePSF = wavePSF/max(wavePSF(:));
        
    % Alpha mask
    mask = ones(size(wavePSF));
    for k = 9:-1:0
        idx = find(wavePSF(:)<(k+1)/100);
        mask(idx) = k/10;
    end
    
    set(app.psfDensityPlotHandle, ...
        'XData', app.components.psf.supportX, ...
        'YData', app.components.psf.supportY, ...
        'CData', wavePSF, ...
        'AlphaData', mask);
    
    for k = 1:size(app.centralConeOutlinesArcMin,1)
        set(app.coneOutlineOnPSFPlotHandles(k), ...
            'XData', app.centralConeOutlinesArcMin(k,1,:), ...
            'YData', app.centralConeOutlinesArcMin(k,2,:));
    end
    
    if (~isempty(app.opticsViewXLimsArcMin))
        set(app.opticsView, ...
            'XLim', app.opticsViewXLimsArcMin, ...
            'YLim', app.opticsViewYLimsArcMin);
    end
    
end



        