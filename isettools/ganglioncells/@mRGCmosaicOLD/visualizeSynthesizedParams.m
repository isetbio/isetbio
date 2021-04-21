function visualizeSynthesizedParams(obj)

    eccDegs = sqrt(sum(obj.synthesizedRFparams.rfCenterPositionDegs.^2,2));
    
    plotParams(201, eccDegs, obj.synthesizedRFparams.visual, 'visual');
    plotParams(202, eccDegs, obj.synthesizedRFparams.retinal, 'retinal'); 
end

function plotParams(figNo, eccDegs, synthParams, domain)
    hFig = figure(figNo);
    set(hFig, 'Name', sprintf('%s params', domain), 'Position', [10 10 700 700]);
    subplot(2,2,1);
    plot(eccDegs, synthParams.surroundCharacteristicRadiiDegs./...
                  synthParams.centerCharacteristicRadiiDegs, 'k.');
    xlabel('ecc (degs)');
    ylabel('surround/center radius');
    set(gca, 'FontSize', 16);
    axis 'square';
    
    subplot(2,2,2);
    plot(eccDegs, synthParams.surroundPeakSensitivities./...
                  synthParams.centerPeakSensitivities, 'k.');
    xlabel('ecc (degs)');
    ylabel('surround/center peak sensitivity');
    set(gca, 'FontSize', 16);
    axis 'square';
    
    subplot(2,2,3);
    integratedCenterSens = synthParams.centerPeakSensitivities .* ...
                           (synthParams.centerCharacteristicRadiiDegs).^2;
    integratedSurroundSens = synthParams.surroundPeakSensitivities .* ...
                           (synthParams.surroundCharacteristicRadiiDegs).^2;
    plot(eccDegs, integratedSurroundSens ./...
                  integratedCenterSens, 'k.');
    xlabel('ecc (degs)');
    ylabel('surround/center integrated sensitivity');
    set(gca, 'FontSize', 16);
    axis 'square';
              
end
