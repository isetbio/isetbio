function str = description(obj, varargin)
% Photopigment text description
%
% Syntax:
%   str = description(obj, varargin)
%
% Description:
%    generate description string for the provided object
%
% Inputs:
%    obj      - The photo pigment objct
%    varargin - (Optional) Any additional information to add? Default empty
%
% Outputs:
%    str      - String containing the photopigment text description.
%
% Optional key/value pairs:
%    None.
%

str = sprintf('Cone aperture:\t %.2f (w) x %.2f (h) um\n', ...
    obj.pdWidth * 1e6, obj.pdHeight * 1e6);
str = addText(str, sprintf('Cone area: \t %.2f um^2\n', obj.area * 1e12));
str = addText(str, sprintf('Gap:\t %.2f(w) x %.2f (h) um\n', ...
    obj.gapWidth * 1e6, obj.gapHeight * 1e6));
str = addText(str, ...
    sprintf('Photopigment density: %.2f (L), %.2f (M), %.2f (S)\n', ...
    obj.opticalDensity(1), obj.opticalDensity(2), obj.opticalDensity(3)));
str = addText(str, ...
    sprintf('Peak efficiency: %.2f (L), %.2f (M), %.2f (S)\n', ...
    obj.peakEfficiency(1), obj.peakEfficiency(2), obj.peakEfficiency(3)));

end