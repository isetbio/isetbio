function [tabulatedDensity, tabulatedEcc] = generateLookUpDensityTables(rfPositions, eccentricitySamplesNum, lambda, neuronalType, whichEye)

    ecc = sqrt(sum(rfPositions .^ 2, 2));
    eccentricitySamplingVector = logspace(log10(lambda), log10(max(ecc)), eccentricitySamplesNum);
    eccentricitySamplingVector  = [-fliplr(eccentricitySamplingVector) 0 eccentricitySamplingVector ];
    [tabulatedEccX, tabulatedEccY] = meshgrid(eccentricitySamplingVector);
    tabulatedEcc = [tabulatedEccX(:) tabulatedEccY(:)];
        
     % Determine smallest spacing (delta)
    switch (neuronalType)
        case 'cone'
            [~,~,tabulatedDensity] = coneSizeReadData(...
                'eccentricity', sqrt(sum(tabulatedEcc.^2,2))*1e-6, ...  % ecc in meters
                'angle', atan2d(tabulatedEccY(:), tabulatedEccX(:)), ...
                'whichEye', whichEye);
        case 'mRGC'
           
        otherwise
            error('Unknown neuronalType: ''%s''.', neuronalType)
    end    
        
    displayTableAs2DMap = false;
    if (displayTableAs2DMap)
        figure()
        X = reshape(tabulatedEcc(:,1), (2*eccentricitySamplesNum+1)*[1 1]);
        Y = reshape(tabulatedEcc(:,2), (2*eccentricitySamplesNum+1)*[1 1]);
        Z = reshape(tabulatedDensity, (2*eccentricitySamplesNum+1)*[1 1]);
        figure
        contourf(X,Y,Z,10)
        colorbar
        axis 'square'
    end
    
end
