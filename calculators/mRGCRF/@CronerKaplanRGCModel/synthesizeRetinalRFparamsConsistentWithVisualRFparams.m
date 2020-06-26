function synthesizedRFParams = synthesizeRetinalRFparamsConsistentWithVisualRFparams(obj, retinalCenterRadiiMicrons, retinalCenterPositionMicrons)

    deconvolutionModel = obj.computeDeconvolutionModel();
    
    surroundToCenterIntegratedSensitivitySigma = 0.075;
    surroundCenterRadiusRatioMean = 0.15;
    surroundCenterRadiusRatioStd = 0.00;
    
    cellsNum = numel(retinalCenterRadiiMicrons);
    
    % Convert RF center positions from retinal microns to visual degs
    retinalEccentricitiesMicrons = sqrt(sum(retinalCenterPositionMicrons.^2,2.0));
    retinalEccentricitiesDegs = (WatsonRGCModel.rhoMMsToDegs(retinalEccentricitiesMicrons/1000.0))';
    
    % Convert RF center radii from retinal microns to visual degs
    pMinus = retinalEccentricitiesMicrons - retinalCenterRadiiMicrons;
    pPlus = retinalEccentricitiesMicrons + retinalCenterRadiiMicrons;
    retinalCenterRadii = (0.5*(WatsonRGCModel.rhoMMsToDegs(pPlus/1000.0) - WatsonRGCModel.rhoMMsToDegs(pMinus/1000.0)))';
    
    % Stats of surround/center radius
    surroundToCenterRadiusRatio = normrnd(surroundCenterRadiusRatioMean, surroundCenterRadiusRatioStd, [1 cellsNum]);
    
    % Preallocate memory
    visualCenterRadii = zeros(1, cellsNum);
    visualSurroundRadii = zeros(1, cellsNum);
    visualCenterPeakSensitivities = zeros(1, cellsNum);
    visualSurroundPeakSensitivities = zeros(1, cellsNum);
    
    retinalSurroundRadii = zeros(1, cellsNum);
    retinalCenterPeakSensitivities = zeros(1, cellsNum);
    retinalSurroundPeakSensitivities = zeros(1, cellsNum);
    
    for cellIndex = 1:cellsNum
        retinalEccDegs = retinalEccentricitiesDegs(cellIndex);
        retinalCenterRadius = retinalCenterRadii(cellIndex);
        [visualCenterRadii(cellIndex),closestEccIndices] = determineVisualRadius(retinalCenterRadius, retinalEccDegs, deconvolutionModel);
        
        % Visual surround radius from visual center radius and surround-to-center radius ratio
        visualSurroundRadii(cellIndex) = visualCenterRadii(cellIndex) / surroundToCenterRadiusRatio(cellIndex);
        
        % Determine the corresponding retinal surround radius
        retinalSurroundRadii(cellIndex) = determineRetinalRadius(visualSurroundRadii(cellIndex), retinalEccDegs, deconvolutionModel, closestEccIndices);
        
        % Compute visual peak sensitivity for center 
        visualCenterPeakSensitivities(cellIndex) = obj.centerPeakSensitivityFunction(obj.centerPeakSensitivityParams, visualCenterRadii(cellIndex));
       
        % Compute visual peak sensitivity for surround based on the stats of integrated sensitivity ratio in the visual data
        surroundToCenterIntegratedSensitivityRatio = normrnd(0.466 + 0.007*retinalEccDegs, surroundToCenterIntegratedSensitivitySigma);
        % Not less than 0.1 and not more than 0.9
        surroundToCenterIntegratedSensitivityRatio = min([max([surroundToCenterIntegratedSensitivityRatio 0.1]) 0.9]);
            
        visualSurroundPeakSensitivities(cellIndex) = surroundToCenterIntegratedSensitivityRatio * visualCenterPeakSensitivities(cellIndex) * (visualCenterRadii(cellIndex)/visualSurroundRadii(cellIndex))^2;
         
        % Determine corresponding retinal gain for center
        gainAttenuation = determineRetinalGainAttenuation(retinalCenterRadii(cellIndex), retinalEccDegs, deconvolutionModel, closestEccIndices);
        retinalCenterPeakSensitivities(cellIndex) = visualCenterPeakSensitivities(cellIndex) / gainAttenuation;
        
        % Determine corresponding retinal gain for surround
        gainAttenuation = determineRetinalGainAttenuation(retinalSurroundRadii(cellIndex), retinalEccDegs, deconvolutionModel, closestEccIndices);
        retinalSurroundPeakSensitivities(cellIndex) = visualSurroundPeakSensitivities(cellIndex) / gainAttenuation;
    end % cellIndex
    
    synthesizedRFParams = struct(...
        'eccDegs', retinalEccentricitiesDegs, ...       % ecc of RGCs within the target patch
        'visual', struct(...                            % visual RF properties (assuming Polans PSFs)
            'centerRadiiDegs', visualCenterRadii, ...
            'surroundRadiiDegs', visualSurroundRadii, ...
            'centerPeakSensitivities', visualCenterPeakSensitivities,...
            'surroundPeakSensitivities', visualSurroundPeakSensitivities...
            ), ...
        'retinal', struct(...                          % retinal RF properties
            'centerRadiiDegs', retinalCenterRadii, ...
            'surroundRadiiDegs', retinalSurroundRadii, ...
            'centerPeakSensitivities', retinalCenterPeakSensitivities,...
            'surroundPeakSensitivities', retinalSurroundPeakSensitivities...
            ) ...
        );
end

function retinalGainAttenuation = determineRetinalGainAttenuation(retinalRadius, retinalEccDegs, deconvolutionModel, closestEccIndices)
    % Compute ecc distances for the closest ecc
    d1 = abs(deconvolutionModel.tabulatedEccentricities(closestEccIndices(1))-retinalEccDegs);
    d2 = abs(deconvolutionModel.tabulatedEccentricities(closestEccIndices(2))-retinalEccDegs);
    
    % Compute attenuations for the 2 closest eccentricities
    gainAttenuation1 =  deconvolutionModel.modelFunctionGain(deconvolutionModel.fittedParamsGain(closestEccIndices(1),:), retinalRadius);
    gainAttenuation2 =  deconvolutionModel.modelFunctionGain(deconvolutionModel.fittedParamsGain(closestEccIndices(2),:), retinalRadius);
    
    % Compute weighted average of retinal gain attenuations
    retinalGainAttenuation = gainAttenuation1 * d2/(d1+d2) + gainAttenuation2 * d1/(d1+d2);
end


function retinalRadius =  determineRetinalRadius(visualRadius, retinalEccDegs, deconvolutionModel, closestEccIndices)
    % Compute ecc distances for the closest ecc
    d1 = abs(deconvolutionModel.tabulatedEccentricities(closestEccIndices(1))-retinalEccDegs);
    d2 = abs(deconvolutionModel.tabulatedEccentricities(closestEccIndices(2))-retinalEccDegs);
    
    % Compute retinal pooling radius functions for the 2 closest eccentricities
    retinalRadiusModel1 = deconvolutionModel.modelFunctionRadius(deconvolutionModel.fittedParamsRadius(closestEccIndices(1),:), visualRadius);
    retinalRadiusModel2 = deconvolutionModel.modelFunctionRadius(deconvolutionModel.fittedParamsRadius(closestEccIndices(2),:), visualRadius);
   
    % Compute weighted average of retinal radii
    retinalRadius = retinalRadiusModel1 * d2/(d1+d2) + retinalRadiusModel2 * d1/(d1+d2);
end


function [visualRadius, closestEccIndices] = determineVisualRadius(retinalRadius, retinalEccDegs, deconvolutionModel)
    % Find closest tabulated eccentricities to this cell's eccentricity
    [~,idx] = sort(abs(abs(deconvolutionModel.tabulatedEccentricities) - retinalEccDegs));
    d1 = abs(deconvolutionModel.tabulatedEccentricities(idx(1))-retinalEccDegs);
    d2 = abs(deconvolutionModel.tabulatedEccentricities(idx(2))-retinalEccDegs);
    closestEccIndices = idx(1:2);
    
    % Compute retinal pooling radius functions for a range of visual
    % radii (up to 10 x retinal radius) for the eccentricity closest to this cell's eccentricity
    visualRadiusSupport = linspace(retinalRadius, retinalRadius*10, 100);
    retinalPoolingRadiusModel1 = deconvolutionModel.modelFunctionRadius(deconvolutionModel.fittedParamsRadius(idx(1),:), visualRadiusSupport);
    retinalPoolingRadiusModel2 = deconvolutionModel.modelFunctionRadius(deconvolutionModel.fittedParamsRadius(idx(2),:), visualRadiusSupport);
        
    % Choose the visual center radius that corresponds to the model
    % retinal pooling radius that is closest to this cells center radius
    [~,index] = min(abs(retinalRadius-retinalPoolingRadiusModel1));
    visualRadius1 = visualRadiusSupport(index);
    [~,index] = min(abs(retinalRadius-retinalPoolingRadiusModel2));
    visualRadius2 = visualRadiusSupport(index);
    visualRadius = visualRadius1 * d2/(d1+d2) + visualRadius2 * d1/(d1+d2);
end

