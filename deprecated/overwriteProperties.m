function target = overwriteProperties(target, source)
% Overwrite fields of a default structure with fields from source
%
%  target = overwriteProperties(target, source)
%
% overwriteProperties runs through all the fields in source, and tries to
% copy them to target. When this fails, we ignore this field.
% 
% 
%
% example:
%
% param.dT = 1;
% defs     = rgcDefaultSimulation('basic');
% defs = overwriteProperties(defs, param);
%

%% Set parameters
warning('This function is deprecated. Use MATLAB command setstructfields');
% errorStr = 'Unrecognized property: %s';

%% Overwrite properties
propertyNames = fieldnames(source);

for ii = 1 : length(propertyNames)
    
    % Check to see if the current field is a legal target
    if isfield(target, propertyNames{ii})
        % Copy field contents to target
        target.(propertyNames{ii}) = source.(propertyNames{ii});
    else
        % Throw error when a failure occurs
        % error(errorStr, propertyNames{ii}); %or not
    end
end

