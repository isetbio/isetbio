function testWatsonModel
    % Unit tests
    RGCmodels.Watson.unitTest.labelValidity();
    RGCmodels.Watson.unitTest.meridianConversions();
    
    % Demo plots
    RGCmodels.Watson.plot.figure1OfWatson2014();
    RGCmodels.Watson.plot.coneDensity2DMaps();
    
    fprintf('All unit tests passed\n');
end
