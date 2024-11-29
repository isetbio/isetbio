% Method to parse the input cal struct
function parseInputCal(obj)
    % load all known fields
    unifiedFieldNames = keys(obj.fieldMap);
    
    % Pass 0: check if fields associated with 'runtimeData__xxx' properties
    % of calStruct class exist in the inputCal, and if they do not add them
    % temporarily with empty values so we can do the parsing/mapping in a
    % unified way.
    if (obj.verbosity > 0)
        fprintf('\n---------------------------------------------------------\n');
        fprintf('<strong>Parsing input cal: Phase 0 (checking for runtime fields) </strong>\n');
    end
    
    % Do not report a warning when accessing S_device, S_ambient
    obj.reportSdev_Samb_warning = false;
    
    for k = 1:numel(unifiedFieldNames)
        % current unified name
        unifiedName = unifiedFieldNames{k};
        
        % retrieve path in input cal
        if (obj.inputCalHasNewStyleFormat)
            calPath = obj.fieldMap(unifiedName).newCalPath;
        else
            calPath = obj.fieldMap(unifiedName).oldCalPath;
        end
        
        % Check if the associated property name is a 'runtimeData' property
        fieldInfo = obj.fieldMap(unifiedName);
        if ~isempty(strfind(fieldInfo.propertyName, 'runtimeData'))
            
            % Check if the runtime field already exists in the inputCal
            addEmptyField = true;
            if (obj.inputCalHasNewStyleFormat)
                dotIndices  = strfind(calPath,'.');
                subStructFieldName = calPath(dotIndices+1:end);
                if ((obj.isFieldOrProperty(obj.inputCal, 'runtimeData')) && ...
                    (obj.isFieldOrProperty(obj.inputCal.runtimeData, subStructFieldName)))
                    addEmptyField = false;
                end
            else
                if obj.isFieldOrProperty(obj.inputCal, calPath)
                    addEmptyField = false;
                end
            end
            
            if (addEmptyField) && (~isempty(calPath))
                % add empty field
                eval(sprintf('obj.inputCal.%s = [];', calPath));
                if (obj.verbosity > 0)
                   fprintf('Added runtime field: %s\n',  calPath);
                end
            end
            
        end
    end
    
    if (obj.verbosity > 0)
        fprintf('<strong>Finished phase 0 parsing. </strong>\n\n');
    end
    
    
    % Pass 1: load properties
    if (obj.verbosity > 0)
        fprintf('\n---------------------------------------------------------\n');
        fprintf('<strong>Parsing input cal: Phase I </strong>\n');
    end
    
    validInputFieldsNum = 0;
    for k = 1:numel(unifiedFieldNames)
        % current unified name
        unifiedName = unifiedFieldNames{k};
        
        % retrieve path in input cal
        if (obj.inputCalHasNewStyleFormat)
            calPath = obj.fieldMap(unifiedName).newCalPath;
        else
            calPath = obj.fieldMap(unifiedName).oldCalPath;
        end
        
        if ~isempty(calPath)
            % get the corresponding private property name
            propertyName = obj.fieldMap(unifiedName).propertyName;
            
            % make sure all sub-structs in the calPath exist in the inputCal
            pathIsValid = true;
            subStruct   = obj.inputCal;
            dotIndices  = strfind(calPath,'.');
            if isempty(dotIndices)
                if ~obj.isFieldOrProperty(subStruct, calPath)
                    if (obj.verbosity > 0)
                        fprintf(2,'>>>> Field   %-25s not found in input cal.\n', sprintf('''%s''',calPath));
                    end
                    %error('>>>> Invalid path for field: ''%s''.\n', calPath);
                    pathIsValid = false;
                end
            else
                p = 1;
                for dotNo = 1:length(dotIndices)
                    subStructFieldName = calPath(p:dotIndices(dotNo)-1);
                    if ~obj.isFieldOrProperty(subStruct, subStructFieldName)
                        if (obj.verbosity > 0)
                            fprintf(2,'>>>> Field   %-25s not found in input cal.\n',  sprintf('''%s''',calPath));
                        end
                        pathIsValid = false;
                        break;
                    end
                    p = dotIndices(dotNo)+1;
                    eval(sprintf('subStruct = subStruct.%s;', subStructFieldName));
                end
                % last field
                subStructFieldName = calPath(p:end);
                if ~obj.isFieldOrProperty(subStruct, subStructFieldName)
                   pathIsValid = false;
                end
            end
            
            if pathIsValid
                propertyValue = eval(sprintf('obj.inputCal.%s;',calPath));
                validInputFieldsNum = validInputFieldsNum + 1;
            else
                propertyValue = [];
                if (obj.verbosity > 1)
                    fprintf(2,'>>>> inputCal does not contain the path ''%s''.\n     Property %s set to [].\n', calPath, propertyName);
                    %eval('inputCalFields = obj.inputCal');
                    fprintf('Hit enter to continue.\n\n');
                    pause;
                end
            end
            
            % and set it
            eval(sprintf('obj.%s = propertyValue;',propertyName));
            
            if (obj.verbosity > 0)
                % Feedback on what hapenned
                if pathIsValid
                    fprintf('%02d. %-45s <- obj.inputCal.%s \n', validInputFieldsNum, propertyName, calPath);
                else
                    fprintf('--- %-45s :: Not found in input cal. Set to [].\n',  propertyName);
                end
            end
        else
            if (obj.verbosity > 1)
                fprintf(2, 'A cal path has not been mapped for unifiedName: ''%s''.\n', unifiedName);
            end
        end
    end
    
    if (obj.verbosity > 0)
        fprintf('<strong>Finished phase I parsing. </strong>\n\n');
        fprintf('<strong>Parsing input cal: Phase II </strong>\n');
    end
    
    % Pass 2: convert any properties that need conversion
    % Note. we have to convert in second pass because fieldnames in a map
    % are not guaranteed to be retrieved in the order they were inserted
    % so if we need a field in the conversion of another field we may not
    % have it.
    
    for k = 1:numel(unifiedFieldNames)
        % current unified name
        unifiedName = unifiedFieldNames{k};
        
        % retrieve path in input cal
        if (obj.inputCalHasNewStyleFormat)
            calPath = obj.fieldMap(unifiedName).newCalPath;
        else
            calPath = obj.fieldMap(unifiedName).oldCalPath;
        end
        
        if ~isempty(calPath)
            % get the corresponding private property name
            propertyName = obj.fieldMap(unifiedName).propertyName;
            
            % check if we need to convert the property to old-style format
            if isfield(obj.fieldMap(unifiedName), 'newToOldConversionFname') && (obj.inputCalHasNewStyleFormat)
                % get conversion function handle
                conversionFunctionHandle = obj.fieldMap(unifiedName).newToOldConversionFname;
                
                if (obj.verbosity > 0)
                    % feedback to the user
                    fprintf('Converting to old-style format: ''obj.inputCal.%s''\n', calPath);
                end
                
                % obtain its current value
                eval(sprintf('propertyValue = obj.%s;',propertyName));
            
                % convert value to old-style format
                propertyValue = conversionFunctionHandle(propertyValue);
                
                % and update the private property
                eval(sprintf('obj.%s = propertyValue;',propertyName));
            end
            
       end
    end
    
    if (obj.verbosity > 0)
        fprintf('<strong>Finished phase II parsing.</strong>\n');
        fprintf('---------------------------------------------------------\n\n');
    end
    
    % Start reporting a warning when accessing S_device, S_ambient
    obj.reportSdev_Samb_warning = true;
    
    
end



