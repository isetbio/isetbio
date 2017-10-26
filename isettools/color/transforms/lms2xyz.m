function xyz = lms2xyz(lms)
% Convert stockman lms to xyz 10 deg
%
% Syntax:
%   xyz = lms2xyz(lms)
%
% Description:
%    Convert stockman lms to xyz 10 deg
%
% Inputs:
%    lms - Probably in RGB image format. *Should check.*
%
% Outputs:
%    xyz - XYZ Formatted image data
%
% Notes:
%    * Why is the Input variable format in question and not already known?
%
% See Also:
%    xyz2lms
%
% (c) ImagEval, 2012

% Examples:
%{
   xyz = lms2xyz(lms)
%}

if notDefined('lms'), error('lms required'); end

if ndims(lms) == 3
    % RGBW format
    xyz = imageLinearTransform(lms, colorTransformMatrix('lms2xyz'));
elseif ismatrix(lms)
    % XW format - Not debugged thoroughly
    if size(lms,1) == 3 && size(lms,2) ~= 3
        xyz = lms' * colorTransformMatrix('lms2xyz');
    elseif size(lms,1) ~= 3 && size(lms,2) == 3
        xyz = lms * colorTransformMatrix('lms2xyz');
    else
        error('Ambiguous lms shape');
    end
end

end
