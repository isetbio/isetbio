% Method to set deepest fields of a struct/substruct by ensuring that
% the fields being set already exist in the struct.
% Usage: params = setExistingFieldsInStruct(params, ...
%               {   'field1', field1Value; ...
%                   'field2', field2Value; ...
%                    ...
%                   'fieldN', 'fieldNValue ...
%                });
% 
%
%   9/15/2016    npc   Wrote it.
%

function struct = setExistingFieldsInStruct(struct, keyValuePairs)
    fieldsNum = size(keyValuePairs,1);
    if (size(keyValuePairs,2) ~= 2)
        error('keyValuePairs must be an Nx2 cell array');
    end
    for k = 1:fieldsNum
        if (isfield(struct, keyValuePairs{k,1}))
            struct.(keyValuePairs{k,1}) = keyValuePairs{k,2};
        else
            error('Struct does not contain a field named: ''%s', keyValuePairs{k,1});
        end
    end
end