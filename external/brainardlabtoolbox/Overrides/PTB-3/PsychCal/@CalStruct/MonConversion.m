% Method to generate the old-style mon parameter  
function mon = MonConversion(obj, propertyValue)
    primariesNum        = size(obj.inputCal.rawData.gammaCurveMeasurements,2); 
    gammaSamples        = size(obj.inputCal.rawData.gammaCurveMeasurements,3); 
    spectralSamples     = size(obj.inputCal.rawData.gammaCurveMeasurements,4); 
    for primaryIndex = 1:primariesNum
        for gammaPointIndex = 1:gammaSamples
            firstSample = (gammaPointIndex-1)*spectralSamples + 1;
            lastSample  = gammaPointIndex*spectralSamples;
            mon(firstSample:lastSample, primaryIndex) = ...
                reshape(squeeze(propertyValue(primaryIndex,gammaPointIndex,:)), ...
                [spectralSamples 1]);
        end
    end
end