% Method to determine if a struct/object contains a field/property with a
% given name
function result = isFieldOrProperty(structOrObject, fieldOrPropertyName)
    if isstruct(structOrObject)
        result = isfield(structOrObject, fieldOrPropertyName);
    elseif isobject(structOrObject)
        result = any(strcmp(properties(structOrObject), fieldOrPropertyName));
    else
        result = false;
    end
end