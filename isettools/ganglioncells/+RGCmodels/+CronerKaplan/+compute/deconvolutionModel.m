function model = deconvolutionModel()
% Assemble the center/surround deconvolution model using data from a deconvolution
% analysis conducted using PSFs at different eccentricities (see
% RGCmodels.CronerKaplan.compute.PolansOpticsDeconvolutionFiles
%

    imposedRefractionErrorDiopters = 0.0;
    meridian = 'horizontal';
    model.tabulatedEccRadiiDegs = sort(...
                             [0.0 0.1 0.2 0.3 0.5 0.7 1.0 ...
                             1.5 2.0 2.5 3.0 4.0 5.0 6.0 ...
                             7.0 8.0 9.0 10  11  12  13 ...
                             14  15  16  17  18  19  20 ...
                             21  22  23  24  25], 'ascend');
                     
    subregions = {'Center', 'Surround'};
    for subregionIndex = 1:2
        subregionName = subregions{subregionIndex};
        
        % Initialize model
        maxConeInputsNumInRFcenter = 30;
        maxSurroundVariations = 20;
        if (strcmpi(subregionName, 'center'))
            model.center.coneInputsNum = nan(numel(model.tabulatedEccRadiiDegs), maxConeInputsNumInRFcenter);
            model.center.characteristicRadiusDegs = nan(numel(model.tabulatedEccRadiiDegs), maxConeInputsNumInRFcenter);
            model.center.peakSensitivity = nan(numel(model.tabulatedEccRadiiDegs), maxConeInputsNumInRFcenter);
        else
            model.surround.centerConeInputsNum = nan(numel(model.tabulatedEccRadiiDegs), maxConeInputsNumInRFcenter);
            model.surround.characteristicRadiusDegs = nan(numel(model.tabulatedEccRadiiDegs), maxConeInputsNumInRFcenter, maxSurroundVariations);
            model.surround.peakSensitivity = nan(numel(model.tabulatedEccRadiiDegs), maxConeInputsNumInRFcenter, maxSurroundVariations);
        end
        
        
        for eccIndex = 1:numel(model.tabulatedEccRadiiDegs)
            % Load data
            eccRadiusDegs = model.tabulatedEccRadiiDegs(eccIndex);
            load(dataFileName(eccRadiusDegs, subregionName, imposedRefractionErrorDiopters, meridian), ...
                'deconvolutionStruct', 'quadrants', 'subjectIDs');
           
            % Get the deconvolution data
            theData = deconvolutionStruct{1,1}.data;
            
            % Parse the deconvolution data for all center cone input configs
            coneInputConfigLabels = keys(theData);

            for coneInputConfigIndex = 1:numel(coneInputConfigLabels)
                
                % Load data for this cone input config
                coneInputConfigLabel = coneInputConfigLabels{coneInputConfigIndex};
                deconvolutionData = theData(coneInputConfigLabel);
                coneInputsNum = str2double(strrep(coneInputConfigLabel, '-coneInput', ''));
                
                % Assemble model
                if (strcmpi(subregionName, 'center'))
                    % Update center model
                    model.center.coneInputsNum(eccIndex,coneInputsNum) = coneInputsNum;
                    model.center.characteristicRadiusDegs(eccIndex,coneInputsNum) = deconvolutionData.characteristicRadiusDegs;
                    model.center.peakSensitivity(eccIndex,coneInputsNum) = deconvolutionData.peakSensitivity;
                else
                    % Update surround model
                    % Retrieve number of computed surround variations
                    surroundVariations = numel(deconvolutionData.nominalSurroundRetinalCharacteristicRadii);
                    
                    model.surround.nominalSurroundRetinalCharacteristicRadii(eccIndex,coneInputsNum,1:surroundVariations) = ...
                        deconvolutionData.nominalSurroundRetinalCharacteristicRadii;
                    model.surround.characteristicRadiusDegs(eccIndex,coneInputsNum,1:surroundVariations) = ...
                        deconvolutionData.characteristicRadiusDegs;
                    model.surround.peakSensitivity(eccIndex,coneInputsNum,1:surroundVariations) = ...
                        deconvolutionData.peakSensitivity;
                end
                
            end % coneInputConfigIndex 
        end % eccIndex
    end % subregionIndex
    
end

function val = dataFileName(eccRadiusDegs, subregion, imposedRefractionErrorDiopters, meridian)

    rootDir = fullfile(RGCmodels.CronerKaplan.constants.centerSurroundDeconvolutionDataDir, sprintf('%smeridian', meridian));
    val = fullfile(rootDir, ...
            sprintf('EccRadius_%2.1fdegs_RefractionError_%2.2fD_%sDeconvolution.mat', ...
            eccRadiusDegs, imposedRefractionErrorDiopters, subregion));
end
