function demoCronerKaplanRGCModel

    ck = CronerKaplanRGCModel();
    
    % synthesis options - use default
    synthesisOptions = ck.synthesisOptions;
%     synthesisOptions = struct( ...
%                 'randomizeCenterRadii', ~true, ...
%                 'randomizeCenterSensitivities', ~true, ...
%                 'randomizeSurroundRadii', ~true, ...
%                 'randomizeSurroundSensitivities', true);
            
    % Simulate CronerKaplan results
    ck.simulateCronerKaplanResults('synthesisOptions', synthesisOptions, ...
        'generateVideo', true, 'generatePlots', true)
    
    

end
