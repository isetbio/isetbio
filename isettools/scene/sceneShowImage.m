function rgb = sceneShowImage(scene, displayFlag, gam)
% Render an image of the scene data
%
% Syntax:
%    rgb = sceneShowImage(scene, [displayFlag], [gam])
%
% Description:
%    The rendering can be either of photons or energy values. This is
%    called from the sceneWindow, so that axes are in that window. If you
%    call this from the command line, a new figure is displayed.
%
% Inputs:
%    scene       - The scene structure
%    displayFlag - (Optional) A boolean value indicating whether or not to
%                  display the scene. Default is 1 (true).
%    gam         - (Optional) The gamma value. Default is 1.
%
%  Outputs:
%    rgb         - The scene data image
%
% Optional key/value pairs:
%    None.
%
% Notes:
%    * TODO: Shouldn't we select the axes for rendering here?  There is
%      only one axis in the scene and oi window. But if we ever go to more,
%      this routine should  say which axis should be used for rendering.
%    * N.B. The source contains executable examples of usage, which can be
%      accessed by typing 'edit sceneShowImage.m' in the command window.
%

% History:
%    xx/xx/03       Copyright ImagEval Consultants, LLC, 2003.
%    12/21/17  jnm  Formatting & fix examples
%    01/26/18  jnm  Formatting update to match Wiki.

% Examples:
%{
    scene = sceneCreate;
    rgb = sceneShowImage(scene);
    sceneShowImage(scene);
    sceneShowImage(scene, 0);
%}

if notDefined('scene'), cla; return;  end
if notDefined('displayFlag'), displayFlag = 1; end
if notDefined('gam'), gam = 1; end

% Force to lower case and no spaces
wList = sceneGet(scene, 'wavelength');
photons = sceneGet(scene, 'photons');
row = sceneGet(scene, 'row'); 
col = sceneGet(scene, 'col');

if isempty(photons)
    cla
    sprintf('ISET Warning:  Data are not available');
    return;
end
    
% This displays the image in the GUI. The displayFlag flag determines how
% imageSPD converts the data into a displayed image. It is set from the
% GUI in the function sceneShowImage.
rgb = imageSPD(photons, wList, gam, row, col, displayFlag);
axis image;
axis off

end