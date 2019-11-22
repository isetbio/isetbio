function [oi, rect] = oiCrop(oi, rect)
% Crop the data in an optical image
%
% Syntax:
%   [oi, rect] = oiCrop(oi, [rect])
%
% Description:
%    Crop the data (photons) in the scene or optical image to within the
%    specified rectangle, rect. If rect is not defined a graphical routine
%    is initiated for selecting the rectangle. The values of rect are also
%    returned.
%
%    Because these are multispectral data, we can't use the usual imcrop.
%    Instead, we use vcROISelect to return the selected data. Then we turn
%    the selected data into a photon image and reset the relevant
%    parameters in the scene structure.
%
%    There are examples contained in the code below. To access, type 'edit
%    oiCrop.m' into the Command Window.
%
% Inputs:
%    oi   - Struct. The optical image structure
%    rect - (Optional) Matrix. The rectangular area to crop
%
% Outputs:
%    oi   - Struct. The modified optical image structure.
%    rect - Matrix. The rectangular area.
%
% Optional key/value pairs:
%    None.
%

% History:
%    xx/xx/03       Copyright ImagEval Consultants, LLC, 2003.
%    03/06/18  jnm  Formatting
%    06/25/19  JNM  Minor formatting changes. Fleshed out example.


% Example:
%{
    % Create a scene
    scene = sceneCreate('uniform');
    scene = sceneSet(scene, 'fov', 15);  % Reasonably large
    scene = sceneAdjustLuminance(scene, 10 ^ -10);

    % Create an optical image
    oi = oiCreate;

    % No lens shading
    optics = oiGet(oi, 'optics');
    optics = opticsSet(optics, 'cos4th', 'off');
    oi = oiSet(oi, 'optics', optics);

    oi = oiCompute(oi, scene);
    oi2 = oiCrop(oi)
%}

if notDefined('oi'), error('You must define an optical image.'); end
if notDefined('rect')
    [roiLocs, rect] = vcROISelect(oi);
else
    cmin = rect(1);
    cmax = rect(1) + rect(3);
    rmin = rect(2);
    rmax = rect(2) + rect(4);
    [c, r] = meshgrid(cmin:cmax, rmin:rmax);
    roiLocs = [r(:), c(:)];
end

wAngular = oiGet(oi, 'wangular');
sz = oiGet(oi, 'size');
% wave = oiGet(oi, 'nwave');

% The number of selected columns and rows
c = rect(3) + 1;
r = rect(4) + 1;

% These are in XW format.
photons = vcGetROIData(oi, roiLocs, 'photons');
photons = XW2RGBFormat(photons, r, c);

oi = oiClearData(oi);
oi = oiSet(oi, 'photons', photons);
newSz = oiGet(oi, 'size');
oi = oiSet(oi, 'wangular', (newSz(2) / sz(2)) * wAngular);

[illuminance, meanIll] = oiCalculateIlluminance(oi);
oi = oiSet(oi, 'illuminance', illuminance);
oi = oiSet(oi, 'meanIlluminance', meanIll);

end