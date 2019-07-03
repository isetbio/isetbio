function fName = oiSaveImage(oi, fName, gam)
% Write PNG image approximating appearance of photon data
%
% Syntax:
%   fName = oiSaveImage(oi, [fName], [gam]);
%
% Description:
%    Save out an RGB image of the photon image as a PNG file. If the name
%    is not passed in, then the user is queried to select the fullpath name
%    of the output file. This routine is used for scenes. sceneSaveImage
%    is used for scenes.
%
% Inputs:
%    oi    - Struct. An optical image structure.
%    fName - (Optional) String. Full file path & name. Default is select an
%            existing rgb image from the current session.
%    gam   - (Optional) Matrix. Matrix containing gamma information.
%
% Outputs:
%    fName - String. Full file name & path of the generated PNG image.
%
% Optional key/value pairs:
%    None.
%

% History:
%    xx/xx/03       Copyright ImagEval Consultants, LLC, 2003.
%    03/06/18  jnm  Formatting
%    06/24/19  JNM  Minor formatting adjustments

if ~exist('oi', 'var') || isempty(oi), oi = vcGetObject('oi'); end

% Get RGB file name (tif)
if notDefined('fName')
    fName = vcSelectDataFile('stayput', 'w', 'png', 'Image file (png)');
end

% Get rgb image from photon data. Gamma either defined here or from the
% open window.
if ~exist('gam', 'var')
    RGB = oiGet(oi, 'rgb image');
else
    RGB = oiGet(oi, 'rgb image', gam);
end

imwrite(RGB, fName, 'png');

end