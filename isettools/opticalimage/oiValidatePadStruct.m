function pad = oiValidatePadStruct(pad)

% Define valid pad values
validPadValues = {'mean photons', 'zero photons'};

% Ensure that the pad is a struct
assert(isstruct(pad), '''pad'' argument must be a struct');


% Ensure that if pad.value exists, it has a valid value
if (isfield(pad, 'value'))
    assert(ischar(pad.value), 'pad.value must be a string');
    assert(ismember(pad.value, {'mean photons', 'zero photons', 'border photons'}), ...
        sprintf('valid values for pad.value are: {''%s'' ''%s''}', validPadValues{1}, validPadValues{2}));
end
end
