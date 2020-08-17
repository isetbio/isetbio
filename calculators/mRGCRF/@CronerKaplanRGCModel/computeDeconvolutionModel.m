function deconvolutionModel = computeDeconvolutionModel(obj, deconvolutionOpticsParams)
    
     
    % Validate the deconvolutionOpticsParams
    obj.validateDeconvolutionOpticsParams(deconvolutionOpticsParams);
    
    defocusMode = 'subjectDefault';
    if (strcmp(defocusMode, 'subjectDefault'))
        imposedRefractionErrorDiopters = 0; 
    else
        imposedRefractionErrorDiopters = 0.01; 
    end
    
     % correction factor for a single cone RF
    x = linspace(-100, 100, 1000);
    coneApertureDiam = 100;
    sigma = coneApertureDiam/2/3;
    y = exp(-(abs(x')/sigma).^2)  * exp(-(abs(x)/sigma).^2);
    y50 = exp(-(abs(x')/sigma).^50)  * exp(-(abs(x)/sigma).^50);
    correctionFactorForDifferenceBetweenFlatopAndGaussianArea = ...
        sum(y50(:))/sum(y(:));
    
    deconvolutionModel.correctionFactorForDifferenceBetweenFlatopAndGaussianArea = correctionFactorForDifferenceBetweenFlatopAndGaussianArea;
    
    
    deconvolutionModel.center = computeCenterDeconvolutionModel(obj,  imposedRefractionErrorDiopters, deconvolutionOpticsParams);
    deconvolutionModel.surround = computeSurroundDeconvolutionModel(obj,  imposedRefractionErrorDiopters, deconvolutionOpticsParams);
    
    
    
 
    
end

function deconvolutionModel = computeCenterDeconvolutionModel(obj,  imposedRefractionErrorDiopters, deconvolutionOpticsParams)
    
    deconvolutionQuadrant = deconvolutionOpticsParams.quadrantsToCompute{1};
    deconvolutionSubject = deconvolutionOpticsParams.PolansWavefrontAberrationSubjectIDsToCompute(1);
    
    deconvolutionModel = computeCenterDeconvolutionModelForSpecificQuadrantAndSubject(obj, ...
        imposedRefractionErrorDiopters, deconvolutionQuadrant, deconvolutionSubject);
end

function deconvolutionModel = computeSurroundDeconvolutionModel(obj,  imposedRefractionErrorDiopters, deconvolutionOpticsParams)
    
    deconvolutionQuadrant = deconvolutionOpticsParams.quadrantsToCompute{1};
    deconvolutionSubject = deconvolutionOpticsParams.PolansWavefrontAberrationSubjectIDsToCompute(1);
    
    deconvolutionModel = computeSurroundDeconvolutionModelForSpecificQuadrantAndSubject(obj, ...
        imposedRefractionErrorDiopters, deconvolutionQuadrant, deconvolutionSubject);
end

function deconvolutionModel = computeSurroundDeconvolutionModelForSpecificQuadrantAndSubject(obj, ...
         imposedRefractionErrorDiopters, deconvolutionQuadrant, deconvolutionSubject)
    
    deconvolutionModel.subjectID = deconvolutionSubject;
    deconvolutionModel.quadrant = deconvolutionQuadrant;
    deconvolutionModel.tabulatedEccentricityRadii = obj.deconvolutionEccs;
    
    for eccIndex = 1:numel(deconvolutionModel.tabulatedEccentricityRadii)

        % Load deconvolution file for this eccentricity
        dataFileName = fullfile(obj.psfDeconvolutionDir,...
            sprintf('ecc_-%2.1f_surroundDeconvolutions_refractionError_%2.2fD.mat', ...
            deconvolutionModel.tabulatedEccentricityRadii(eccIndex), imposedRefractionErrorDiopters));
        load(dataFileName, 'deconvolutionStruct', 'quadrants', 'subjectIDs');
        
         assert((numel(subjectIDs) == 1) && (subjectIDs == deconvolutionSubject), ...
            sprintf('Deconvolution file does not contain subject %d', deconvolutionSubject));
        
        assert((numel(quadrants) == 1) && (strcmp(quadrants{1},deconvolutionQuadrant)), ...
            sprintf('Deconvolution file does not contain quadrant''%s''', deconvolutionQuadrant));
        
        % Load the data for the RFsurround
        theData = deconvolutionStruct{1,1}.surround.data;
        theDataLabels = keys(theData);
        
        visualGains = [];
        retinalGains = [];
        visualCharacteristicRadius = [];
        retinalCharacteristicRadius = [];
        for centerConeInputConfigIndex = 1:numel(theDataLabels)
            centerConeInputConfig = theDataLabels{centerConeInputConfigIndex};
            deconvolutionData = theData(centerConeInputConfig);

           % deconvolutionModel.surroundRetinalNominalCharacteristicRadius(eccIndex, coneInputsNum, :) = deconvolutionData.nominalSurroundRetinalCharacteristicRadii;
            visualGains(centerConeInputConfigIndex,:) = deconvolutionData.visualGains;
            retinalGains(centerConeInputConfigIndex,:) = deconvolutionData.retinalGains;
            visualCharacteristicRadius(centerConeInputConfigIndex,:) = deconvolutionData.visualSigmas;
            retinalCharacteristicRadius(centerConeInputConfigIndex,:) = deconvolutionData.retinalSigmas;
        end % coneInputConfigIndex
        

        deconvolutionModel.visualGain(eccIndex,:,:) = visualGains;
        deconvolutionModel.retinalGain(eccIndex,:,:) = retinalGains;
        deconvolutionModel.visualCharacteristicRadius(eccIndex,:,:)  = visualCharacteristicRadius;
        deconvolutionModel.retinalCharacteristicRadius(eccIndex,:,:) = retinalCharacteristicRadius;
        
    end % eccIndex
    
end


function deconvolutionModel = computeCenterDeconvolutionModelForSpecificQuadrantAndSubject(obj, ...
        imposedRefractionErrorDiopters, deconvolutionQuadrant, deconvolutionSubject)
    
    deconvolutionModel.subjectID = deconvolutionSubject;
    deconvolutionModel.quadrant = deconvolutionQuadrant;
    deconvolutionModel.tabulatedEccentricityRadii = obj.deconvolutionEccs;
    
      
    for eccIndex = 1:numel(deconvolutionModel.tabulatedEccentricityRadii)

        % Load deconvolution file for this eccentricity
        dataFileName = fullfile(obj.psfDeconvolutionDir,...
            sprintf('ecc_-%2.1f_centerDeconvolutions_refractionError_%2.2fD.mat', ...
            deconvolutionModel.tabulatedEccentricityRadii(eccIndex), imposedRefractionErrorDiopters));
        load(dataFileName, 'deconvolutionStruct', 'quadrants', 'subjectIDs');
        
        assert((numel(subjectIDs) == 1) && (subjectIDs == deconvolutionSubject), ...
            sprintf('Deconvolution file does not contain subject %d', deconvolutionSubject));
        
        assert((numel(quadrants) == 1) && (strcmp(quadrants{1},deconvolutionQuadrant)), ...
            sprintf('Deconvolution file does not contain quadrant''%s''', deconvolutionQuadrant));
        
        % Load the data for the RFcenter
        theData = deconvolutionStruct{1,1}.center.data;
        theDataLabels = keys(theData);
        
        for coneInputConfigIndex = 1:numel(theDataLabels)
            coneInputConfig = theDataLabels{coneInputConfigIndex};
            deconvolutionData = theData(coneInputConfig);
            coneInputsNum = str2double(strrep(coneInputConfig, '-coneInput', ''));

            deconvolutionModel.centerConeInputsNum(eccIndex,coneInputsNum) = coneInputsNum;
            deconvolutionModel.visualGain(eccIndex,coneInputsNum) = deconvolutionData.visualGain;
            deconvolutionModel.retinalGain(eccIndex,coneInputsNum) = deconvolutionData.retinalGain;
            deconvolutionModel.visualCharacteristicRadiusMin(eccIndex,coneInputsNum) = deconvolutionData.minVisualSigma;
            deconvolutionModel.visualCharacteristicRadiusMax(eccIndex,coneInputsNum) = deconvolutionData.maxVisualSigma;
            deconvolutionModel.retinalCharacteristicRadiusMin(eccIndex,coneInputsNum) = deconvolutionData.minRetinalSigma;
            deconvolutionModel.retinalCharacteristicRadiusMax(eccIndex,coneInputsNum) = deconvolutionData.maxRetinalSigma;
            deconvolutionModel.retinalCharacteristicRadiusMax(eccIndex,coneInputsNum) = deconvolutionData.maxRetinalSigma;
        end
    end
end

% 
% 
% 
% function deconvolutionModel = OLDcomputeSurroundDeconvolutionModel(obj, tabulatedEccentricities, imposedRefractionErrorDiopters, deconvolutionOpticsParams)
%     subjectsToAverage = deconvolutionOpticsParams.PolansWavefrontAberrationSubjectIDsToAverage;
%     quadrantsToAverage = deconvolutionOpticsParams.quadrantsToAverage;
%     
%     % Use WatsonRGCModel to retrieve cone apertures along the nasal meridian
%     % for the tabulated eccentricities
%     w = WatsonRGCModel();
%     coneRFSpacingsDegs  = w.coneRFSpacingAndDensityAlongMeridian(abs(tabulatedEccentricities), ...
%             'nasal meridian','deg', 'deg^2', ...
%             'correctForMismatchInFovealConeDensityBetweenWatsonAndISETBio', false);
%     % Cone aperture is a percentage of the cone spacing
%     coneApertureRadii = WatsonRGCModel.coneApertureToDiameterRatio * 0.5 * coneRFSpacingsDegs;
%     coneApertureRadiusAtZeroDegs = coneApertureRadii(1);
%     
%     for eccIndex = 1:numel(tabulatedEccentricities)
%         % Eccentricity
%         eccDegs = tabulatedEccentricities(eccIndex);
%         
%         % Load deconvolution file for this eccentricity
%         dataFileName = fullfile(obj.psfDeconvolutionDir,...
%             sprintf('ecc_-%2.1f_deconvolutions_refractionError_%2.2fD.mat', eccDegs, imposedRefractionErrorDiopters));
%         load(dataFileName, 'retinalPoolingRadii', 'visualRadius', 'visualGain', 'subjectIDs', 'quadrants');
%         retinalPoolingRadiiOriginal = retinalPoolingRadii;
%         
%         % Compute the minimum possible RF center radius as the mean (over all subjects)
%         % convolution of the cone aperture at 0 deg with the PSF at 0 deg 
%         if (eccDegs == 0)
%             % Find the 2 closest retinal pooling radius
%             [dd,idx] = sort(abs( coneApertureRadiusAtZeroDegs-retinalPoolingRadii));
%             % Compute the minimal visual radius as the weighed average of the
%             % visual radii corresponding to the 2 closest retinal pooling radii
%             d1 = dd(2)/(dd(1)+dd(2));
%             d2 = dd(1)/(dd(1)+dd(2));
%             visualRadii = d1 * squeeze(visualRadius(idx(1),:,:)) + d2 * squeeze(visualRadius(idx(2),:,:));
%             minVisualRadiusDegs = median(median(visualRadii));
%         end
% 
%         % Get data for the quadrant of interest
%         visualRadius = CronerKaplanRGCModel.quadrantData(visualRadius, quadrantsToAverage, quadrants, subjectsToAverage, subjectIDs);
%         visualGain = CronerKaplanRGCModel.quadrantData(visualGain, quadrantsToAverage, quadrants, subjectsToAverage, subjectIDs);
%         
% 
%         % We will fit the relationship between visual and retinal radii
%         % using a saturating function, but only for retinal radii >= cone aperture radii
%         % So select these data here:
%         idx = find(retinalPoolingRadii >= coneApertureRadii(eccIndex));
%         retinalPoolingRadii = retinalPoolingRadii(1,idx);
%         visualRadius = visualRadius(:,idx);
%         visualGain = visualGain(:,idx);
%         
%         % Compute the median visual radius over subjects/ecc quadrants
%         medianVisualRadius = (median(visualRadius,1, 'omitnan'))';
%         
%         % Compute the median gain over all subjects/ecc quadrants
%         medianVisualGain = (median(visualGain,1, 'omitnan'))';
%        
%         % Fit the relationship: retinalRadius(visualRadius) = Model(visualRadius, params)
%         [modelFunctionRadius, fittedParamsRadius(eccIndex,:)] = fitRadiusData(medianVisualRadius, retinalPoolingRadii);
%         % Fit the relationship: visualGain(retinalRadius) = Model(retinalRadius, params)
%         [modelFunctionGain, fittedParamsGain(eccIndex,:)] = fitGainData(medianVisualGain, retinalPoolingRadii);
%         
%         nonAveragedVisualRadius{eccIndex} = visualRadius;
%         nonAveragedVisualGain{eccIndex} = visualGain;
%     end % eccIndex
%     
%     deconvolutionModel = struct(...
%         'retinalPoolingRadii', retinalPoolingRadiiOriginal, ...
%         'tabulatedEccentricities', tabulatedEccentricities, ...
%         'fittedParamsRadius', fittedParamsRadius, ...
%         'fittedParamsGain', fittedParamsGain, ...
%         'modelFunctionRadius', modelFunctionRadius, ...
%         'modelFunctionGain', modelFunctionGain, ...
%         'minVisualRadiusDegs', minVisualRadiusDegs ...
%         );
%      % Other meta-parameters
%      deconvolutionModel.coneApertureRadii = coneApertureRadii;
%      deconvolutionModel.opticsParams = deconvolutionOpticsParams;
%      deconvolutionModel.nonAveragedVisualRadius = nonAveragedVisualRadius;
%      deconvolutionModel.nonAveragedVisualGain = nonAveragedVisualGain;
% end
% 
% function [modelFunction, fittedParams] = fitRadiusData(visualRadius, retinalRadius)
%     % For better fitting extend the visual radius data to 1.5
%     visual  = [visualRadius; [0.7 0.8 1 1.5]'];
%     retinal = [retinalRadius [0.7 0.8 1 1.5] ]';
%         
%     modelFunction = @(p,x)(p(1) - (p(2)-p(1))*exp(-p(3)*x));
%     initialParams = [5 10 0.2];
%     [fittedParams, fittedParamsSE] = nonLinearFitData(visual, retinal, modelFunction, initialParams);
% end
% 
% function [modelFunction, fittedParams] = fitGainData(visualGain, retinalRadius)
%     % For better fitting extend the visual radius data to 1.5
%     visual  = [visualGain; [1 1 1]'];
%     retinal = [retinalRadius [0.8 1 2]]';
%         
%     modelFunction = @(p,x)((p(1)*x.^p(3))./(x.^p(3)+p(2)));
%     initialParams = [0.9 0.002 2];
%     [fittedParams, fittedParamsSE] = nonLinearFitData(retinal, visual, modelFunction, initialParams);
% end
% 
% function [fittedParams, fittedParamsSE] = nonLinearFitData(x,y, modelFunction, initialParams)
%     opts.RobustWgtFun = 'talwar';
%     opts.MaxIter = 1000;
%     [fittedParams,~,~,varCovarianceMatrix,~] = nlinfit(x,y,modelFunction,initialParams,opts);
%     % standard error of the mean
%     fittedParamsSE = sqrt(diag(varCovarianceMatrix));
%     fittedParamsSE = fittedParamsSE';
% end