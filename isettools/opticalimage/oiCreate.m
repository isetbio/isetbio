function oi = oiCreate(oiType,varargin)
%Create an optical image structure
%
%   oi = oiCreate(oiType,varargin)
%
% The optical image represents the spectral irradiance at the sensor. The
% irradiance is computed from the scene radiance, using the information in
% the optics structure that is attached to this oi.
%
% The oi spectrum is normally inherited from the scene.  To specify a
% spectrum for the optical image use
%      oi = oiCreate('default');
%      oi = initDefaultSpectrum('hyperspectral');
%
% Types of OI structures
%  {'human','default'}  - Human shift-invariant optics based on Marimont
%                   and Wandell (1994, JOSA)
%  {'wvf human'}    - Human shift-invariant optics based on mean
%                   wavefront abberration from Thibos et al. (2009,
%                   Ophthalmic & Physiological Optics)
%  {'diffraction'}  - Diffraction limited optics,f/4, no diffuser or data
%  {'uniformd65'}   - Turns off offaxis to make uniform D65 image
%  {'uniformee'}    - Turns off offaxis to make uniform equal energy image
%
% Example:
%   oi = oiCreate('human');
%   oi = oiCreate('uniform d65');  % Used for lux-sec vs. snr measurements.
%   oi = oiCreate('uniform EE');   %
%   oi = oiCreate('diffraction');
%   oi = oiCreate('wvf human');
%   pupilMM = 4; oi = oiCreate('wvf human',pupilMM)
%
% See also:  sceneCreate, opticsCreate
%
% Copyright ImagEval Consultants, LLC, 2003.

if notDefined('oiType'),  oiType = 'human'; end

% Default is to use the diffraction limited calculation
oi.type = 'opticalimage';
oi.name = vcNewObjectName('opticalimage');  % Get a fresh name
oi = oiSet(oi, 'bit depth', 32);            % Single precision.  Perhaps no longer needed

oiType = ieParamFormat(oiType);
switch oiType 
    case {'default','human', 'mwhuman'}
        % Marimont and Wandell optics, which is a simple shift-invariant
        % but wavelength-dependent model.  This is a little faster than the
        % wvf human, so we made it the default.  They differ a little.
        %
        % oi = oiCreate('human');
        oi = oiCreate('diffraction limited');
        oi = oiSet(oi, 'diffuser method','skip');
        oi = oiSet(oi, 'consistency',1);
        oi = oiSet(oi, 'optics', opticsCreate('human'));
        oi = oiSet(oi, 'name','human-MW');
        oi = oiSet(oi, 'lens', Lens('wave', oiGet(oi, 'wave')));
        
    case {'wvfhuman','shiftinvariant'}
        % A human lens specified using the WVF toolbox method
        % oi = oiCreate('wvf human',pupilMM,zCoefs,wave)
        
        oi = oiCreate('diffraction limited');
        oi = oiSet(oi, 'diffuser method','skip');
        oi = oiSet(oi, 'consistency',1);
        oi = oiSet(oi, 'optics', opticsCreate('wvf human',varargin{:}));
        oi = oiSet(oi, 'name','human-WVF');
        oi = oiSet(oi, 'lens', Lens('wave', oiGet(oi, 'wave')));
        
    case {'diffractionlimited','diffraction'}
        % Default optics is f# = 4, diffraction limited
        optics = opticsCreate('diffraction limited');
        oi = oiSet(oi,'optics',optics);
        oi = oiSet(oi, 'name','diffraction');

        % Set up the default glass diffuser with a 2 micron blur circle,
        % but skipped
        oi = oiSet(oi, 'diffuser method','skip');
        oi = oiSet(oi, 'diffuser blur', 2*10^-6);
        oi = oiSet(oi, 'consistency', 1);
        
    case {'uniformd65'}
        % Uniform, D65 optical image.  No cos4th falloff, huge field of
        % view (120 deg). Used in lux-sec SNR testing and scripting
        oi = oiCreateUniformD65;
               
    case {'uniformee'}
        % Uniform, equal energy optical image.  No cos4th falloff, huge
        % field of view (120 deg). 
        oi = oiCreateUniformEE;
        
    otherwise
        error('Unknown oiType');
end

end

%--------------------------------------------
function oi = oiCreateUniformD65
%
%  Create a uniform, D65 image with a very large field of view.  The
%  optical image is created without any cos4th fall off so it can be used
%  for lux-sec SNR testing.
%

% This does not yet extend in the IR, but it should.  See notes in
% sceneCreate.
% This does not yet extend in the IR, but it should.  See notes in
% sceneCreate.
scene = sceneCreate('uniform d65');
scene = sceneSet(scene,'hfov',120);

oi = oiCreate('diffraction');
oi = oiSet(oi,'optics fnumber',1e-3);   % Basically perfect optics
oi = oiCompute(scene,oi);
% vcAddObject(oi); oiWindow;

end

%--------------------------------------------
function oi = oiCreateUniformEE
%
%  Create a uniform, equal energy image with a very large field of view.
%  The optical image is created without any cos4th fall off so it can be
%  used for lux-sec SNR testing.
%

% This does not yet extend in the IR, but it should.  See notes in
% sceneCreate.
scene = sceneCreate('uniform ee');
scene = sceneSet(scene,'hfov',120);

oi = oiCreate('diffraction');
oi = oiSet(oi,'optics fnumber',1e-3);   % Basically perfect optics
oi = oiCompute(scene,oi);
% vcAddObject(oi); oiWindow;

end
