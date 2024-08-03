% Method to determine whether the inputCal has the expected basic
% fields.
function determineInputCalValidity(obj)
    % an input cal must have at least a 'describe', and a 'rawdata'
    % section. Check for these and exit if they do not exist.
    
    if (~isfield(obj.inputCal, 'describe'))
        error('inputCal does not contain a ''describe'' part.\n');
    end
        
    if (obj.inputCalHasNewStyleFormat)
        if (~isfield(obj.inputCal, 'rawData'))
            error('inputCal does not contain a ''rawData'' part.\n');
        end
    else
        if (~isfield(obj.inputCal, 'rawdata'))
            error('inputCal does not contain a ''rawdata'' part.\n');
        end
    end
end

        