function allViews(app)
    CSFGeneratorApp.render.coneMosaicView(app, 'initialize');
    CSFGeneratorApp.render.roiView(app, 'initialize');
    CSFGeneratorApp.render.opticsView(app, 'initialize');
    CSFGeneratorApp.render.stimulusView(app, 'initialize');
    CSFGeneratorApp.render.psychometricFunctionView(app, 'initialize');
    CSFGeneratorApp.render.psychometricFunctionView(app, 'update');
    CSFGeneratorApp.render.csfView(app, 'initialize');
    CSFGeneratorApp.render.csfView(app, 'update');
end
