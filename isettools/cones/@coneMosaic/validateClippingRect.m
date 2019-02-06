function validateClippingRect(clippingRect)
% Validate the clippingRect
%
% Syntax:
%   validateClippingRect(clippingRect)
%
% Description
%    The clippingRect must be a structure with fields 'xo', yo', 'width', 
%    and 'height', all defined in visual degrees.
%
% Inputs:
%    clippingRect - The clipping rect, a struct(see above)
%
% Outputs:
%    None
%
% Optional key/value pairs:
%    None.

% History:
%    2/6/19  NPC  ISETBIO Team, 2019


    validFnamesString = '''xo'', ''yo'', ''width'' and ''height'', all specified in degrees';
    validFnames = {'xo', 'yo', 'width', 'height'};
    if (~isstruct(clippingRect))
        error('ClippingRect must be a struct with fields: %s', validFnamesString);
    end
    
    fNames = fieldnames(clippingRect);
    specifiedFieldNames = {};
    for k = 1:numel(fNames)
        if (~ismember(fNames{k},  validFnames))
            error('Fieldname ''%s'' is not valid. Valid field names: %s', fNames{k}, validFnamesString);
        else
            specifiedFieldNames{numel(specifiedFieldNames)+1} = fNames{k};
        end
    end
    
    for k = 1:numel(validFnames)
        if (~ismember(validFnames{k}, specifiedFieldNames))
            error('Fieldname ''%s'' was not specified. Clipping rect must contain %s', validFnames{k}, validFnamesString);
        end
    end
end