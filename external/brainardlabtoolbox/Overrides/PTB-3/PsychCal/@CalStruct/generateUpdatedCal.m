
% Method to generate a cal struct (of the old format)
% by parsing the obj.oldFormatFieldMap.
function cal = generateUpdatedCal(obj)

    if (obj.verbosity > 0)
        fprintf('Generating cal with old-style format.\n');
    end
    
    cal = struct();
    
    % Get all the mapped unified field names
    unifiedFieldNames = keys(obj.fieldMap);
    
    % Do not generate warnings when S_device, S_ambient are accessed
    obj.reportSdev_Samb_warning = false;
    
    for k = 1:numel(unifiedFieldNames)
        calPath = obj.fieldMap(unifiedFieldNames{k}).oldCalPath;
        if ~isempty(calPath)
            propertyName  = obj.fieldMap(unifiedFieldNames{k}).propertyName;
            propertyValue = eval(sprintf('obj.%s;',propertyName));
            if (obj.verbosity > 0)
                fprintf('%02d. Loading cal.%-30s <-- %s\n',k, calPath, propertyName); 
            end
            eval(sprintf('cal.%s = propertyValue;',calPath));
        end
    end
    
    % Start generating warnings when S_device, S_ambient are accessed
    obj.reportSdev_Samb_warning = true;
    
    % Mark the cal as exported from NewStyleCalStruct
    cal.describe.isExportedFromNewStyleCalStruct = true;  
end
