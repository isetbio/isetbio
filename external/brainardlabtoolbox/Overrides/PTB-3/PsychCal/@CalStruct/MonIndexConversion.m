% Method to generate the old-style monIndex parameter  
function monIndex = MonIndexConversion(obj, propertyValue)
    trialsNum           = size(obj.inputCal.rawData.gammaCurveMeasurements,1);
    primariesNum        = size(obj.inputCal.rawData.gammaCurveMeasurements,2); 
    for trialIndex = 1:trialsNum 
        for primaryIndex = 1:primariesNum   
            monIndex{trialIndex, primaryIndex} = squeeze(propertyValue(trialIndex, primaryIndex,:));  
        end
    end
end