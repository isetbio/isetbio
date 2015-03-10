function xyz = lms2xyz(lms)
% Convert stockman lms to xyz 10 deg
%
%    xyz = lms2xyz(lms)
%
% Input:
%   lms:  Probablyk in RGB image format.  Should check.
% Example:
%
% See also:  xyz2lms
%
% (c) ImagEval, 2012

if notDefined('lms'), error('lms required'); end

if ndims(lms) == 3
    % RGBW format
    xyz = imageLinearTransform(lms, colorTransformMatrix('lms2xyz'));
elseif ismatrix(lms)
    % XW format - Not debugged thoroughly
    if (size(lms,1) == 3) && size(lms,2) ~= 3
        xyz = lms * colorTransformMatrix('lms2xyz');
    elseif    (size(lms,1) ~= 3) && size(lms,2) == 3
        xyz = lms' * colorTransformMatrix('lms2xyz');
    else
        error('Ambiguous lms shape');
    end
end

end

