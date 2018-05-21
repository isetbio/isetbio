function rgb = displayShowImage(d, displayFlag, varargin)
% Render an image of the display image data
%
% Syntax:
%   rgb = displayShowImage(d, [displayFlag], [varargin])
%
% Description:
%    Render an image of the display image data.
%
% Inputs:
%    d           - Struct. An isetbio display structure.
%    displayFlag - (Optional) Boolean. Boolean indicating whether or not
%                  to display the image. Default is 1(true).
%     varargin   - (Optional) Additional arguments as desired. The first
%                  would be an axes object to use in image.
%
% Outputs:
%    rgb         - Matrix. RGB Image data.
%
% Optional key/value pairs:
%    None.
%

% History:
%    05/xx/14  HJ   Created May, 2014
%    05/15/18  jnm  Formatting

% Examples:
%{
    % ETTBSkip - skipping since 'ha' is not defined.
    d = displayCreate('LCD-Apple');
    rgb = displayShowImage(d);
    displayShowImage(d);
    rgb = displayShowImage(d, 1, ha);
%}
%{
    d = displayCreate('LCD-Apple');
    rgb = displayShowImage(d);
    displayShowImage(d);
    rgb = displayShowImage(d, 1);
%}

if isempty(d), cla; return; end
if notDefined('displayFlag'), displayFlag = 1; end
if ~isempty(varargin), ha = varargin{1}; end

% Get parameters
wave = displayGet(d, 'wave');
global vcSESSION;

% Main display image
rgb = displayGet(d, 'main image');
if isempty(rgb)
    % HJ did this originally. But it should go away.
    if isfield(vcSESSION, 'imgData')
        rgb = vcSESSION.imgData;
    else
        cla;
        warning('No image data found');
        return;
    end
end

% Compute photon image
scene = sceneFromFile(rgb, 'rgb', [], d);
img = sceneGet(scene, 'photons');

% This displays the image in the GUI. The displayFlag flag determines
% how imageSPD converts the data into a displayed image. It is set from
% the GUI in the function displayShowImage.
if ~notDefined('ha'), axes(ha); end

gam = 1;
rgb = imageSPD(img, wave, gam, [], [], displayFlag);
axis image;
axis off

return;
