% Method to check the validity of the requested unified field name.
function isValid = fieldNameIsValid(obj, unifiedFieldName)
%
    isValid = true;
    % make sure fieldname is in the FieldMap
    if (~ismember(unifiedFieldName, keys(obj.fieldMap)))
        isValid = false;
    end
end