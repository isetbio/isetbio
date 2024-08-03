% Getter method for a given unified field name
function fieldValue = get(obj, unifiedFieldName)
	fieldValue = NaN;
    if (~ischar(unifiedFieldName))
        error('>>> Field name must be entered in quotes.\n');
    end
    if (obj.fieldNameIsValid(unifiedFieldName))
        % Find the corresponding property name
        propertyName = obj.fieldMap(unifiedFieldName).propertyName;
        % Call the getter for that property
        fieldValue = eval(sprintf('obj.%s;',propertyName));
    else
        fprintf(2, '>>> Unknown unified field name (''%s''). Cannot get its value.\n', unifiedFieldName);
        obj.printMappedFieldNames(); 
        error('Read denied from unknown field.\n');
    end     
end 