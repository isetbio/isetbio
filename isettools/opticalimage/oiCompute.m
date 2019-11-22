function oi = oiCompute(scene, oi, opticsModel)
% Gateway routine for optical image irradiance calculation
%
% Syntax:
%   oi = oiCompute(oi, scene, [opticsModel])
%
% Description:
%    The spectral irradiance image, on the sensor plane, just before sensor
%    capture is the optical image. This spectral irradiance distribution
%    depends on the scene and the the optics attached to the optical image
%    structure, oi.
%
%    Three types of optical calculations are implemented in ISET. These are
%    selected from the interface in the oiWindow, or they can be set via
%    optics = opticsSet(optics, 'model', parameter)
%
%    This oiCompute function calls the relevant model depending on the
%    setting in opticsGet(optics, 'model');
%
%    The first model is based on diffraction-limited optics and calculated
%    in opticsDLCompute. This computation depends on the diffraction
%    limited parameters (f/#) but little else.
%
%    The second model is shift-invariant optics. This depends on the user
%    setting up a wavelength-dependent OTF defined and included in the
%    optics structure.
%
%    You may specify the shift-invariant optics model to be 'human'. In
%    that case, this program calls the routine humanOI, which uses the
%    human OTF calculations in a shift-invariant fashion.
%
%    You can also build shift invariant optics models specifying wavefront
%    aberrations using Zernike polynomials and the wavefront (wvf<TAB>)
%    methods.
%
%    N.B. You can modify the waitbars with the following commands:
%        Use ieSessionSet('waitbar', 0) to turn off waitbar displays
%        Use ieSessionSet('waitbar', 1) to turn on waitbar displays
%
%    There are examples contained in the code below. To access, type 'edit
%    oiCompute.m' into the Command Window.
%
% Inputs:
%    scene       - Struct. A scene structure.
%    oi          - Struct. An optical image structure.
%    opticsModel - (Optional) String. An optics model. Default is to
%                  retrieve from the provided optical image. Options
%                  include: 'diffractionlimited' and 'shiftInvariant'.
%
% Outputs:
%    oi          - Struct. The modified optical image structure, containing
%                  the calculated irradiance.
%
% Optional key/value pairs:
%    None.
%
% See Also:
%   opticsGet, opticsSet, wvfCreate, opticsCreate
%

% History:
%    xx/xx/05       Copyright ImagEval Consultants, LLC, 2005
%    03/06/18  jnm  Formatting
%    04/07/18  dhb  Leave the gun, keep the canoli. (This means,
%                   deleted broken example, kept the working one.)
%    06/25/19  JNM  Minor formatting adjustments

% Examples:
%{
    scene = sceneCreate;
    oi = oiCreate('human');
    oi = oiCompute(scene, oi);
%}

if notDefined('scene'), error('Scene required.'); end
if notDefined('oi'), error('Optical image required.'); end
if strcmp(oi.type, 'scene') && strcmp(scene.type, 'opticalimage')
    tmp = scene;
    scene = oi;
    oi = tmp;
    clear tmp
end
if notDefined('opticsModel'), opticsModel = oiGet(oi, 'optics model'); end
oi = oiSet(oi, 'distance', sceneGet(scene, 'distance'));

% Compute according to the selected model
switch ieParamFormat(opticsModel)
    case {'diffractionlimited', 'dlmtf'}
        % The skip case is handled by the DL case
        oi = opticsDLCompute(scene, oi);
    case 'shiftinvariant'
        oi = opticsSICompute(scene, oi);
    otherwise
        error('Unknown optics model')
end

% Indicate scene it is derived from
oi = oiSet(oi, 'name', sceneGet(scene, 'name'));

% After computing, the oi data representation should have the same wave as
% the scene
oi = oiSet(oi, 'wave', sceneGet(scene, 'wave'));

% Pad the scene depth map and attach it to the oi.  The padded values are
% set to 0, though perhaps we should pad them with the mean distance.
oi = oiSet(oi, 'depth map', oiPadDepthMap(scene));

end