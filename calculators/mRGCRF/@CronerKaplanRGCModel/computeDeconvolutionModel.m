function deconvolutionModel = computeDeconvolutionModel(obj, deconvolutionOpticsParams)
    
    % Validate the deconvolutionOpticsParams
    obj.validateDeconvolutionOpticsParams(deconvolutionOpticsParams);
    
    defocusMode = 'subjectDefault';
    if (strcmp(defocusMode, 'subjectDefault'))
        imposedRefractionErrorDiopters = 0; 
    else
        imposedRefractionErrorDiopters = 0.01; 
    end

 
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
    
    maxConeInputsNumInRFcenter = 30;
    maxSurroundVariations = 20;
    deconvolutionModel.centerConeInputsNum = nan(numel(deconvolutionModel.tabulatedEccentricityRadii),maxConeInputsNumInRFcenter);
    deconvolutionModel.nominalSurroundRetinalCharacteristicRadii = nan(numel(deconvolutionModel.tabulatedEccentricityRadii),maxConeInputsNumInRFcenter,maxSurroundVariations);
    deconvolutionModel.characteristicRadiusDegs = nan(numel(deconvolutionModel.tabulatedEccentricityRadii),maxConeInputsNumInRFcenter,maxSurroundVariations);
    deconvolutionModel.peakSensitivity = nan(numel(deconvolutionModel.tabulatedEccentricityRadii),maxConeInputsNumInRFcenter,maxSurroundVariations);
      
        
    for eccIndex = 1:numel(deconvolutionModel.tabulatedEccentricityRadii)

        % Load deconvolution file for this eccentricity
        dataFileName = obj.deconvolutionDataFileName(...
            -deconvolutionModel.tabulatedEccentricityRadii(eccIndex), imposedRefractionErrorDiopters, 'surround');
        load(dataFileName, 'deconvolutionStruct', 'quadrants', 'subjectIDs');
        
        assert((numel(subjectIDs) == 1) && (subjectIDs == deconvolutionSubject), ...
            sprintf('Deconvolution file does not contain subject %d', deconvolutionSubject));
        
        assert((numel(quadrants) == 1) && (strcmp(quadrants{1},deconvolutionQuadrant)), ...
            sprintf('Deconvolution file does not contain quadrant''%s''', deconvolutionQuadrant));
        
        % Load the data for the RFsurround
        theData = deconvolutionStruct{1,1}.data;
        theDataLabels = keys(theData);
        
        for centerConeInputConfigIndex = 1:numel(theDataLabels)
            coneInputConfig = theDataLabels{centerConeInputConfigIndex};
            deconvolutionData = theData(coneInputConfig);
            coneInputsNum = str2double(strrep(coneInputConfig, '-coneInput', ''));
            
            deconvolutionModel.centerConeInputsNum(eccIndex,coneInputsNum) = coneInputsNum;
            surroundVariationsForThisCenterConeInputConfig = numel(deconvolutionData.nominalSurroundRetinalCharacteristicRadii);
            deconvolutionModel.nominalSurroundRetinalCharacteristicRadii(eccIndex,coneInputsNum,1:surroundVariationsForThisCenterConeInputConfig) = ...
                deconvolutionData.nominalSurroundRetinalCharacteristicRadii;
            deconvolutionModel.characteristicRadiusDegs(eccIndex,coneInputsNum,1:surroundVariationsForThisCenterConeInputConfig) = ...
                deconvolutionData.characteristicRadiusDegs;
            deconvolutionModel.peakSensitivity(eccIndex,coneInputsNum,1:surroundVariationsForThisCenterConeInputConfig) = ...
                deconvolutionData.peakSensitivity;
        end % coneInputConfigIndex
        
    end % eccIndex
    
end


function deconvolutionModel = computeCenterDeconvolutionModelForSpecificQuadrantAndSubject(obj, ...
        imposedRefractionErrorDiopters, deconvolutionQuadrant, deconvolutionSubject)
    
    deconvolutionModel.subjectID = deconvolutionSubject;
    deconvolutionModel.quadrant = deconvolutionQuadrant;
    deconvolutionModel.tabulatedEccentricityRadii = obj.deconvolutionEccs;
    
    maxConeInputsNumInRFcenter = 30;
    deconvolutionModel.centerConeInputsNum = nan(numel(deconvolutionModel.tabulatedEccentricityRadii),maxConeInputsNumInRFcenter);
    deconvolutionModel.characteristicRadiusDegs = nan(numel(deconvolutionModel.tabulatedEccentricityRadii),maxConeInputsNumInRFcenter);
    deconvolutionModel.peakSensitivity = nan(numel(deconvolutionModel.tabulatedEccentricityRadii),maxConeInputsNumInRFcenter);
        
    for eccIndex = 1:numel(deconvolutionModel.tabulatedEccentricityRadii)
        % Load deconvolution file for this eccentricity
        dataFileName = obj.deconvolutionDataFileName(...
            -deconvolutionModel.tabulatedEccentricityRadii(eccIndex), imposedRefractionErrorDiopters, 'center');
        load(dataFileName, 'deconvolutionStruct', 'quadrants', 'subjectIDs');
        
        assert((numel(subjectIDs) == 1) && (subjectIDs == deconvolutionSubject), ...
            sprintf('Deconvolution file ''%s'' does not contain subject %d', dataFileName, deconvolutionSubject));
        
        assert((numel(quadrants) == 1) && (strcmp(quadrants{1},deconvolutionQuadrant)), ...
            sprintf('Deconvolution file ''%s'' does not contain quadrant''%s''', dataFileName, deconvolutionQuadrant));
        
        % Load the data for the RFcenter
        deconvolutionModel.metaData{eccIndex} = deconvolutionStruct{1,1}.metaData;
        theData = deconvolutionStruct{1,1}.data;
        theDataLabels = keys(theData);

        for coneInputConfigIndex = 1:numel(theDataLabels)
            coneInputConfig = theDataLabels{coneInputConfigIndex};
            deconvolutionData = theData(coneInputConfig);
            coneInputsNum = str2double(strrep(coneInputConfig, '-coneInput', ''));

            deconvolutionModel.centerConeInputsNum(eccIndex,coneInputsNum) = coneInputsNum;
            deconvolutionModel.characteristicRadiusDegs(eccIndex,coneInputsNum) = deconvolutionData.characteristicRadiusDegs;
            deconvolutionModel.peakSensitivity(eccIndex,coneInputsNum) = deconvolutionData.peakSensitivity;
            
        end
    end
end