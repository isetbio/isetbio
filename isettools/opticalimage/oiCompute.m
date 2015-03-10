function oi = oiCompute(scene,oi,opticsModel)
% Gateway routine for optical image irradiance calculation
%
%    oi = oiCompute(scene,oi,[opticsModel])
%
% The spectral irradiance image, on the sensor plane, just before sensor
% capture is the optical image.  This spectral irradiance distribution
% depends on the scene and the the optics attached to the optical image
% structure, oi. 
%
% Three types of optical calculations are implemented in ISET.  These are
% selected from the interface in the oiWindow, or they can be set via
% optics = opticsSet(optics,'model', parameter)
%
% This oiCompute function calls the relevant model depending on the setting
% in opticsGet(optics,'model'); 
%
% The first model is based on diffraction-limited optics and calculated in
% opticsDLCompute. This computation depends on the diffraction limited
% parameters (f/#) but little else.
%
% The second model is shift-invariant optics.  This depends on having a
% wavelength-dependent OTF defined and included in the optics structure. 
%
% The third model is a full ray trace model that allows shift-varying 
% wavelength-dependent point spread functions.  This depends on having the
% ray trace information included in the optics model.  The geometric
% information, relative illumination, and blurring calculations are usually
% imported from a full ray trace program, such as Zemax.
%
% It is possible to specify the model in the argument. If it is not
% specified, then the optics data in the oi determines the model.
%
% Finally, there is a special case for computing the human optical image.
% You may specify the optics model to be 'human'.  In that case, this
% program calls the routine humanOI, which uses the human OTF calculations
% in a shift-invariant fashion.
%
% See also: opticsGet for a description of the optics models.
%
% Use ieSessionSet('waitbar',0) to turn off waitbar displays
% Use ieSessionSet('waitbar',1) to turn on waitbar displays
%
% Example
%   oi = oiCompute(scene,oi);
%
%   oi = vcGetObject('oi'); scene = vcGetObject('scene');
%   load siZemaxExample00
%   optics = opticsSet(optics,'model','shiftinvariant');
%   oi = oiSet(oi,'optics',optics);
%   oi = oiCompute(scene,oi);
% 
% Programming note:  I should have written this as oiCompute(oi,scene, ..)
% I now allow this ordering (rather than oiCompute(scene,oi, ...)) by a
% type check at the front that will flip the arguments.
%
% Copyright ImagEval Consultants, LLC, 2005

if notDefined('scene'), error('Scene required.'); end
if notDefined('oi'), error('Optical image required.'); end
if strcmp(oi.type,'scene') && strcmp(scene.type,'opticalimage')
    % disp('oiCompute: flipping arguments')
    tmp = scene; scene = oi; oi = tmp; clear tmp
end
if notDefined('opticsModel'), 
    optics = oiGet(oi,'optics'); 
    opticsModel = opticsGet(optics,'model'); 
end
oi = oiSet(oi, 'distance', sceneGet(scene, 'distance'));

% Compute according to the selected model
opticsModel = ieParamFormat(opticsModel);
switch lower(opticsModel)
    case {'diffractionlimited','dlmtf','skip'}
        % The skip case is handled by the DL case
        oi = opticsDLCompute(scene,oi);
    case 'shiftinvariant'
        oi = opticsSICompute(scene,oi);
    case 'raytrace'
        oi = opticsRayTrace(scene,oi);
    otherwise
        error('Unknown optics model')
end

% Indicate scene it is derived from
oi = oiSet(oi,'name',sceneGet(scene,'name'));

% Pad the scene dpeth map and attach it to the oi.   The padded values are
% set to 0, though perhaps we should pad them with the mean distance.
oi = oiSet(oi,'depth map',oiPadDepthMap(scene));

end