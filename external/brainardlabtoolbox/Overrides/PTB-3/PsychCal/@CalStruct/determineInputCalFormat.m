% Method to determine whether the inputCal has new-style format.
function determineInputCalFormat(obj)
%
    if (~isfield(obj.inputCal,'describe'))
    	error('Input calibration structure has no ''describe'' field.');
    end

    obj.inputCalHasNewStyleFormat = false;
    if isfield(obj.inputCal.describe, 'driver')
        if strcmp(obj.inputCal.describe.driver, 'object-oriented calibration')
            if (~((isfield(obj.inputCal.describe, 'isExportedFromNewStyleCalStruct')) && (obj.inputCal.describe.isExportedFromNewStyleCalStruct == true)))
                obj.inputCalHasNewStyleFormat = true;
            end
        end
    end
end