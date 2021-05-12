function coneMosaicView(app, mode)
    switch (mode)
        case 'initialize'
            initializeConeMosaicView(app);
        case 'update'
            updateConeMosaicViewWithNewData(app);
    end
end

function initializeConeMosaicView(app)

    set(app.coneMosaicView, 'LineWidth', 0.5);
    box(app.coneMosaicView, 'off');
    
    % Do not show the interactions toolbax
    app.coneMosaicView.Toolbar.Visible = 'off';
            
    % Only allow zooming
    app.coneMosaicView.Interactions = [panInteraction zoomInteraction];
            
    % Add listener to zoom-events
    addlistener(app.coneMosaicView, {'XLim', 'YLim'}, 'PostSet', @app.handleConeMosaicViewZoomEvent);
end

function updateConeMosaicViewWithNewData(app)

    if (~isempty(app.coneMosaicViewXLimsDegs))
        switch (app.coneMosaicParams.visualizationDomain)
           case 'degrees'
                domainVisualizationLimits(1:2) = app.coneMosaicViewXLimsDegs;
                domainVisualizationLimits(3:4) = app.coneMosaicViewYLimsDegs;
           case 'microns'
                domainVisualizationLimits(1:2) = app.coneMosaicViewXLimsMicrons;
                domainVisualizationLimits(3:4) = app.coneMosaicViewXLimsMicrons;
        end
    else
        domainVisualizationLimits = [];
    end
    
    mosaicVisualizationView = 'retinal view';
    
    if (~isempty(app.components.coneMosaic.coneRFpositionsDegs))
        app.components.coneMosaic.visualize(...
            'figureHandle', app.mainView, ...
            'axesHandle', app.coneMosaicView, ...
            'visualizationView', mosaicVisualizationView, ...
            'labelCones', true, ...
            'crossHairsOnFovea', true, ...
            'visualizedConeAperture', 'geometricArea', ...
            'domain', app.coneMosaicParams.visualizationDomain, ...
            'domainVisualizationLimits', domainVisualizationLimits, ...
            'noYLabel', true, ...
            'fontSize', 14, ...
            'backgroundColor', 'none', ... %[0.3 0.3 0.3], ...
            'plotTitle',  ' ');
    else
        cla(app.coneMosaicView)
    end
    
    
end
