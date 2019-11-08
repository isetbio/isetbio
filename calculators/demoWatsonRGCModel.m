function demoWatsonRGCModel
    WatsonRGCCalc = WatsonRGCModel('generateAllFigures', true);
    
    
    fprintf('peak cone density: %3.2fk cones/mm2\n', WatsonRGCCalc.peakConeDensity('cones per mm2')/1000);
    fprintf('peak cone density: %3.2fk cones/deg2\n', WatsonRGCCalc.peakConeDensity('cones per deg2')/1000);
    
    
end

