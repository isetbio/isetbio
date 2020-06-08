function synthesizeData(obj,eccDegs,synthesisOptions)
    
    % Noise-free parameters
    centerRadii = obj.centerRadiusFunction(obj.centerRadiusParams, eccDegs);
    surroundRadii = obj.surroundRadiusFunction(obj.surroundRadiusParams, eccDegs);
    centerSurroundRadiusRatios = centerRadii ./ surroundRadii;
    
    centerPeakSensitivities = obj.centerPeakSensitivityFunction(obj.centerPeakSensitivityParams, centerRadii);
    surroundPeakSensitivities = obj.surroundPeakSensitivityFunction(obj.surroundPeakSensitivityParams, surroundRadii);
    surroundCenterPeakSensitivityRatios = surroundPeakSensitivities./centerPeakSensitivities;
    

    maxAttemptsNo = 1000;
    pW = 0.15;
    
    % Synthesize params for each rfUnit
    parfor rfUnit = 1:numel(eccDegs)
        
        meanIntegratedSurroundToCenterSensitivityRatio = 0.466 + eccDegs(rfUnit)*0.007*0.5;
        integratedSurroundToCenterSensitivityRatio = Inf;
        attemptsNo = 0;
        p = abs(randn*pW);
        centerRadii(rfUnit) = 0;
        centerPeakSensitivities(rfUnit) = 0;
        surroundRadii(rfUnit) = 0;
        surroundPeakSensitivities(rfUnit) = 0;
        
        % Adjust surround radii/sensitivity so as to bring the integrated
        % surround/center sensitivity within the distribution of Figure 11
        while (attemptsNo < maxAttemptsNo) && (...
              ((abs(integratedSurroundToCenterSensitivityRatio-meanIntegratedSurroundToCenterSensitivityRatio) > p) || ...
              (integratedSurroundToCenterSensitivityRatio<0.1) || ...
              (integratedSurroundToCenterSensitivityRatio>0.9) || ...
              (centerRadii(rfUnit) < 0) || (centerPeakSensitivities(rfUnit)<0) || ...
              (surroundRadii(rfUnit)< 0) || (surroundPeakSensitivities(rfUnit) < 0)))
              
            
            [centerRadii(rfUnit), centerPeakSensitivities(rfUnit), ...
                surroundRadii(rfUnit), surroundPeakSensitivities(rfUnit)] = drawSurroundAgain(obj, ...
                        eccDegs(rfUnit), centerSurroundRadiusRatios(rfUnit), surroundCenterPeakSensitivityRatios(rfUnit), synthesisOptions);
            
            integratedSurroundToCenterSensitivityRatio = ...
                (surroundRadii(rfUnit)/centerRadii(rfUnit))^2 * (surroundPeakSensitivities(rfUnit)/centerPeakSensitivities(rfUnit));
            
            p = abs(randn*pW);
            attemptsNo = attemptsNo + 1;
            
        end
        
        if (attemptsNo == maxAttemptsNo)
            fprintf('Failed to meet integrated sensitivity ratio after %d attempts for eccentricity %f\n', maxAttemptsNo, eccDegs(rfUnit));
        end
    end % parfor
    
    assert(all(centerPeakSensitivities>0), 'Found center peak sensitivities < 0');
    assert(all(surroundPeakSensitivities>0), 'Found surround peak sensitivities < 0');
    assert(all(surroundRadii>0), 'Found surround radii < 0');
    assert(all(centerRadii>0), 'Found center radii < 0');
    
    obj.synthesizedData = struct(...
        'eccDegs', eccDegs, ...
        'centerRadii', centerRadii, ...
        'surroundRadii', surroundRadii, ...
        'centerPeakSensitivities', centerPeakSensitivities, ...
        'surroundPeakSensitivities', surroundPeakSensitivities ...
        );
end

     
function [centerRadius,  centerPeakSensitivity, surroundRadius, surroundPeakSensitivity] = ...
        drawSurroundAgain(obj, eccDegs, centerSurroundRadiusRatio, surroundCenterPeakSensitivityRatio, synthesisOptions)
    
     if (synthesisOptions.randomizeSurroundRadii)
        % Derive from model of surround size
        surroundRadiusParamsNoisy = normrnd(obj.surroundRadiusParams, obj.surroundRadiusParamsSE);
        surroundRadius = obj.surroundRadiusFunction(surroundRadiusParamsNoisy, eccDegs);
        % Derive from center radius and stats of surround/center radius ratio
        %stochasticRatio = normrnd(obj.centerData('size').radiusRatioToSurroundStats(1), obj.centerData('size').radiusRatioToSurroundStats(2));
        %surroundRadius = centerRadius / stochasticRatio;
     else
        surroundRadius = obj.surroundRadiusFunction(obj.surroundRadiusParams, eccDegs); 
     end
    
     if (synthesisOptions.randomizeCenterRadii)
        centerRadiusParamsNoisy = normrnd(obj.centerRadiusParams, obj.centerRadiusParamsSE);
        centerRadius = obj.centerRadiusFunction(centerRadiusParamsNoisy, eccDegs);
     else
        %centerRadius = surroundRadius * centerSurroundRadiusRatio;
        centerRadius = obj.centerRadiusFunction(obj.centerRadiusParams, eccDegs);
     end

     if (~synthesisOptions.randomizeSurroundRadii) && (synthesisOptions.randomizeCenterRadii)
          surroundRadius = centerRadius / centerSurroundRadiusRatio;
     end
    
    if (synthesisOptions.randomizeCenterSensitivities)
        % Noisy center peak sensitivities
        centerPeakSensitivityParamsNoisy = normrnd(obj.centerPeakSensitivityParams, obj.centerPeakSensitivityParamsSE);
        centerPeakSensitivity = obj.centerPeakSensitivityFunction(centerPeakSensitivityParamsNoisy, centerRadius);
    else
        centerPeakSensitivity = obj.centerPeakSensitivityFunction(obj.centerPeakSensitivityParams, centerRadius);
    end
            
    

    if (synthesisOptions.randomizeSurroundSensitivities)
        surroundPeakSensitivityParamsNoisy = normrnd(obj.surroundPeakSensitivityParams, obj.surroundPeakSensitivityParamsSE);
        surroundPeakSensitivity = obj.centerPeakSensitivityFunction(surroundPeakSensitivityParamsNoisy, surroundRadius);
    else
        % surround sensitivities from noisy center sensitivities based on noise-free ratios
        surroundPeakSensitivity = centerPeakSensitivity * surroundCenterPeakSensitivityRatio;
    end
    
end
