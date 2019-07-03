function pad = oiValidatePadStruct(pad)
% Validate the pad structure
%
% Syntax:
%   pad = oiValidatePadStruct(pad)
%
% Description:
%    Validate and return the existing pad structure.
%
% Inputs:
%    pad - Struct. A structure with a key named value, and a corresponding
%          string value of either 'mean photons' or 'zero photons'.
%
% Outpits:
%    pad - Struct. The validated pad structure.
%
% Optional key/value pairs:
%    None.
%

% Define valid pad values
validPadValues = {'mean photons', 'zero photons'};

% Ensure that the pad is a struct
assert(isstruct(pad), '''pad'' argument must be a struct');

% Ensure that if pad.value exists, it has a valid value
if (isfield(pad, 'value'))
    assert(ischar(pad.value), 'pad.value must be a string');
    assert(ismember(pad.value, {'mean photons', 'zero photons'}), ...
        sprintf('valid values for pad.value are: {''%s'' ''%s''}', ...
        validPadValues{1}, validPadValues{2}));
end

end
