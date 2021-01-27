function XYZ = lms2xyz(lms)
% Convert stockman lms to xyz 10 deg
%
% Syntax:
%   XYZ = lms2xyz(lms)
%
% Description:
%    Convert Stockman lms to XYZ 10 degree
%
%    This function contains examples of usage inline. To access these, type
%    'edit lms2xyz.m' in the Command Window.
%
% Inputs:
%    lms - Matrix. The Stockman LMS, in either RGB or XW format.
%
% Outputs:
%    XYZ - Matrix. The XYZ in RGB format.
%
% Optional key/value pairs:
%    None.
%
% Notes:
%    * [Note - DHB: Someday might be nice to have units arg that specified
%      10 versus 2 degree, perhaps separately for Stockman-Sharpe and XYZ,
%      and perhaps using key/value pairs.]
%
% See Also:
%   xyz2lms
%

% History:
%    xx/xx/12       (c) ImagEval
%    11/01/17  jnm  Comments, formatting, and adding actual example.
%    11/02/17  dhb  This tried to be clever about accepting non XW format.
%                   Took that out.
%                   Extended example to test XW format and remove notes
%                   about usefulness of doing so.
%    11/17/17  jnm  Formatting
%    07/16/19  JNM  Formatting update

% Examples:
%{
    scene = sceneCreate('reflectance chart');
    vcAddAndSelectObject(scene);
    sceneWindow
    imgLMS = sceneGet(scene, 'lms');
    [xwLMS, m, n] = RGB2XWFormat(imgLMS);
    imgXYZ = lms2xyz(imgLMS);
    xwXYZ = lms2xyz(xwLMS);
    imgXYZCheck = XW2RGBFormat(xwXYZ, m, n);
    if max(abs(imgXYZ(:) - imgXYZCheck(:)))
        fprintf('Oops. Something wrong with format conversion\n');
    end
    vcNewGraphWin;
    imagescRGB(lms2srgb(imgXYZ));
%}

if notDefined('lms'), error('lms required'); end

if ndims(lms) == 3
    % RGBW format
    XYZ = imageLinearTransform(lms, colorTransformMatrix('lms2xyz'));
elseif ismatrix(lms)
    if size(lms, 2) ~= 3
        error(['Passed lms not in RGB or XW format. ' ...
               'You might need to transpose input for proper XW']);
    else
        XYZ = lms * colorTransformMatrix('lms2xyz');
    end
end

end
