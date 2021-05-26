function stimulusView(app, mode)

    switch (mode)
        case 'initialize'
            initializeStimulusView(app);
        case 'update'
            updateStimulusViewWithNewData(app);
    end
end

function initializeStimulusView(app)
    app.stimulusPlotHandle = image(app.stimulusView, [-0.01 0.01], [-0.01 0.01], [1 1; 1 1]);
    app.stimulusTextHandle = text(app.stimulusView, 0, 0, '', 'FontSize', 14);
    set(app.stimulusView, 'XTick', [], 'YTick', [], 'XColor', 'none', 'YColor', 'none', 'FontSize', 1);
    axis(app.stimulusView, 'image');
    
    % Do not show the interactions toolbax
    app.stimulusView.Toolbar.Visible = 'off';
            
    % No interactions
    app.stimulusView.Interactions = [];
   
end

function updateStimulusViewWithNewData(app)
    frames = numel(app.products.demoStimulusSceneSequence);

    for iFrame = 1:numel(frames)
        theCurrentScene = app.products.demoStimulusSceneSequence{iFrame};
        rgbImage = sceneGet(theCurrentScene, 'rgb');
        xAxis = 1:size(rgbImage,2);
        yAxis = 1:size(rgbImage,1);
        set(app.stimulusPlotHandle, ...
            'XData', xAxis, ...
            'YData', yAxis, ...
            'CData', rgbImage);
        set(app.stimulusTextHandle, 'Position', [xAxis(1+round(0.02*numel(xAxis))), yAxis(1+round(0.03*numel(yAxis)))]);
        set(app.stimulusTextHandle, 'String', sprintf('frame %d/%d', iFrame, numel(frames)), 'FontSize', 14);
        axis(app.stimulusView, 'image');
        drawnow;
    end
    
end
