% Method to generate the old-style background-dependence spectra parameter
function bgSpectra = BGspectraConversion(obj, propertyValue)
    spectralSamples       = size(obj.inputCal.rawData.gammaCurveMeasurements,4);
    backgroundSettingsNum = size(obj.inputCal.backgroundDependenceSetup.bgSettings,2);
    targetSettingsNum     = size(obj.inputCal.backgroundDependenceSetup.settings,2);
    for backgroundSettingsIndex = 1:backgroundSettingsNum
        tmp = zeros(spectralSamples,targetSettingsNum); 
        
        if (isempty(propertyValue))
            tmp = [];
        else
            for targetSettingsIndex = 1: targetSettingsNum
                tmp(:, targetSettingsIndex) = ...
                reshape(squeeze(propertyValue(backgroundSettingsIndex, targetSettingsIndex, :)), ...
                [spectralSamples  1] );
            end
        end
        
        bgSpectra{backgroundSettingsIndex} = tmp;
    end 
end