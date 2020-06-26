function [tabulatedDensity, tabulatedSpacing, tabulatedEcc] = generateLookUpDensityTables(rfPositions, eccentricitySamplesNum, lambda, neuronalType, whichEye)

    ecc = sqrt(sum(rfPositions .^ 2, 2));
    eccentricitySamplingVector = logspace(log10(lambda), log10(max(ecc)), eccentricitySamplesNum);
    eccentricitySamplingVector  = [-fliplr(eccentricitySamplingVector) 0 eccentricitySamplingVector ];
    [tabulatedEccX, tabulatedEccY] = meshgrid(eccentricitySamplingVector);
    tabulatedEcc = [tabulatedEccX(:) tabulatedEccY(:)];
        
     % Determine smallest spacing (delta)
    switch (neuronalType)
        case 'cone'
            tabulatedDensity = coneDensityFunctionFull(tabulatedEcc, whichEye);
        case 'mRGC'
            tabulatedDensity = mRGCRFDensityFunctionFull(tabulatedEcc, whichEye);
        otherwise
            error('Unknown neuronalType: ''%s''.', neuronalType)
    end    
        
    % Spacing (in microns) from density
    tabulatedSpacing = 1e3 * WatsonRGCModel.spacingFromDensity(tabulatedDensity);
    
    displayTableAs2DMap = ~true;
    if (displayTableAs2DMap)
        figure(123); clf;
        subplot(1,2,1)
        X = reshape(tabulatedEcc(:,1), (2*eccentricitySamplesNum+1)*[1 1]);
        Y = reshape(tabulatedEcc(:,2), (2*eccentricitySamplesNum+1)*[1 1]);
        Z = reshape(tabulatedDensity, (2*eccentricitySamplesNum+1)*[1 1]);
        contourf(X,Y,Z,10)
        title('density');
        colorbar
        axis 'square'
        
        subplot(1,2,2)
        X = reshape(tabulatedEcc(:,1), (2*eccentricitySamplesNum+1)*[1 1]);
        Y = reshape(tabulatedEcc(:,2), (2*eccentricitySamplesNum+1)*[1 1]);
        Z = reshape(tabulatedSpacing, (2*eccentricitySamplesNum+1)*[1 1]);
        contourf(X,Y,Z,10)
        title('spacing');
        colorbar
        axis 'square'
        drawnow;
    end
    
end
