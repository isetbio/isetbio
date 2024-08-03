% Method to make old-style meter type. In OOC calibration we store
% the meter model as a string.
function meterID = MeterTypeConversion(obj, propertyValue)
    % In the new-style cal we store the meter model as a string,
    % whereas in the old-style cal we stored a number only.
    meterID = 0;
    if (strcmp(obj.describe___meterModel, 'PR-650'))
        meterID = 1;
    elseif (strcmp(obj.describe___meterModel, 'PR-655'))
        meterID = 4;
    elseif (strcmp(obj.describe___meterModel, 'PR-670'))
        meterID = 5;
    end
end