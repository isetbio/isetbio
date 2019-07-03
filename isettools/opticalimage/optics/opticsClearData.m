function obj = opticsClearData(obj)
% Clear OTF data from optics. Used after parameter change to optics or OI.
%
% Syntax:
%   obj = opticsClearData(obj)
%
% Description:
%    Clear the data from the optics fields after one of the optics
%    parameters has been changed.
%
%    The input object can be either an optical image or an optics object.
%
% Inputs:
%    obj - Object. The input object.
%
% Outputs:
%    obj - Object. The modified object.
%
% Optional key/value pairs:
%    None.
%

% History:
%    xx/xx/03       Copyright ImagEval Consultants, LLC, 2003.
%    03/09/18  jnm  Formatting
%    06/27/19  JNM  Formatting update

switch lower(obj.type)
    case {'opticalimage', 'oi'}
        optics = oiGet(obj, 'optics');

        if checkfields(optics, 'OTF', 'OTF')
            optics = opticsSet(optics, 'OTFdata', []);
        end
        obj = oiSet(obj, 'optics', optics);

    case 'optics'
        obj = opticsSet(obj, 'OTFdata', []);

    otherwise
        error('Unknown object type.');
end

end