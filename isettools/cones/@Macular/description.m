function str = description(obj, varargin)
% Generate description string for this macular pigment
%
% Syntax:
%   str = description(obj, varargin)
%
% Description:
%    generate description string for the provided object
%
% Inputs:
%    obj      - The macular pigment objct
%    varargin - (Optional) Any additional information to add? Default empty
%
% Outputs:
%    str      - String containing the macular pigment text description.
%
% Optional key/value pairs:
%    None.
%

str = sprintf('Macular pigment density: %.2f\n', obj.density);

end