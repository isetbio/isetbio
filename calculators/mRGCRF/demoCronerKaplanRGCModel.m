function demoCronerKaplanRGCModel

    ck = CronerKaplanRGCModel();
    
    
    % synthesis options - use default
    synthesisOptions = ck.synthesisOptions;
    synthesisOptions = struct( ...
                'randomizeCenterRadii', ~true, ...
                'randomizeCenterSensitivities', ~true, ...
                'randomizeSurroundRadii', ~true, ...
                'randomizeSurroundSensitivities', true);
            
    % Simulate CronerKaplan results
    ck.simulateCronerKaplanResults('synthesisOptions', synthesisOptions, ...
        'generateVideo', true, 'generatePlots', true)
    pause
    
    % Ecc
    eccDegs = logspace(log10(0.01), log10(50), 2000);
    % default syntehsis options

    ck.synthesizeData(eccDegs, synthesisOptions);
    [hFigParams, hFigRatios] = ck.plotSynthesizedData();
    
    ck.plotlabOBJ.exportFig(hFigParams, 'pdf', sprintf('synthesizedParams'), pwd());
    ck.plotlabOBJ.exportFig(hFigRatios, 'pdf', sprintf('synthesizedRatios'), pwd());
        
end
