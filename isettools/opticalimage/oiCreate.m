function [oi,val] = oiCreate(oiType,val,optics,addObject,varargin)
%Create an optical image structure.
%
%   [oi,val] = oiCreate(oiType,val,optics,addObject,varargin)
%
% If val is passed in, the optical image is set to be the number val.
% Otherwise, a new number is selected.
%
% If optics is passed in, this is  attached to the optical image. Otherwise
% the default optics (diffraction limited) are used.
%
% By default, the new optical image is not added to the set of optical
% image objects stored in vcSESSION.  If you want it added, addObject = 1.
% Normally we add it when ready with a vcAddAndSelectObject
%
% OI types include: default, uniformD65, uniformEE.  The latter two are
% used only for lux-snr testing and related.  Almost always we simply
% create a default optical image with a diffraction-limited lens attached.
%
% The spectrum is not set in this call because it is normally inherited
% from the scene.  To specify a spectrum for the optical image use
%      oi = oiCreate('default');
%      oi = initDefaultSpectrum('hyperspectral');
%
% Types of OI:
%  {'default'}    - Diffraction limited optics,f/4, no diffuser or data
%  {'uniformd65'} - Turns off offaxis to make uniform D65 image
%  {'uniformee'}  - Turns off offaxis and creates uniform equal energy image
%
%  {'human'}      - Human shift-invariant optics based on Marimont
%                   and Wandell (1994, JOSA)
%  {'wvf human'}  - Human shift-invariant optics based on mean wavefront
%                   abberration from Thibos et al. (2009, Ophthalmic & Physiological
%                   Optics)
% Example:
%   oi = oiCreate('default');
%   oi = oiCreate('uniform d65');          % Used for lux-sec vs. snr measurements.
%   oi = oiCreate('uniform EE',[],[],0);   % Create an object but don't put it in the vcSESSION
%   oi = oiCreate('uniformEE',[],[],0,(380:4:1068));
%   oi = oiCreate('human');
%   oi = oiCreate('wvf human');
%
% See also:  sceneCreate
%
% Copyright ImagEval Consultants, LLC, 2003.

if notDefined('oiType'),  oiType = 'human'; end
if notDefined('val'),     val = vcNewObjectValue('OPTICALIMAGE'); end
if notDefined('optics'),  optics = opticsCreate('default'); end

% We used to automatically add created OI objects to the list.  Stopped
% doing this July, 2012
if notDefined('addObject'), addObject = 0; end

% Default is to use the diffraction limited calculation
oi.type = 'opticalimage';
oi.name = vcNewObjectName('opticalimage');
oi = oiSet(oi, 'bit depth', 32);  % Single precision.

oiType = ieParamFormat(oiType);
switch oiType 
    case {'default','human', 'mwhuman'}
        % Marimont and Wandell optics, which is a simple shift-invariant
        % but wavelength-dependent model.  This is a little faster than the
        % wvf human, so we made it the default.  They differ a little.
        %
        % oi = oiCreate('human');
        oi = oiCreate('diffraction limited');
        oi = oiSet(oi, 'diffuserMethod','skip');
        oi = oiSet(oi, 'consistency',1);
        oi = oiSet(oi, 'optics', opticsCreate('human'));
        oi = oiSet(oi, 'name','human-MW');
        
    case {'wvfhuman','shiftinvariant'}
        % A human lens specified using the WVF toolbox method
        % oi = oiCreate('wvf human',pupilMM,zCoefs)
        % This is in the class of shift invariant models
        
        % Defaults and then adjust for varargin
        wave = 400:10:700; wave = wave(:);
        pupilMM = 3;
        zCoefs = wvfLoadThibosVirtualEyes(pupilMM);
        if ~isempty(varargin), pupilMM = varargin{1}; end
        if length(varargin) > 1, zCoefs = varargin{2}; end
        if length(varargin) > 2, wave = varargin{3}; end
        
        % Create wavefront parameters
        wvfP = wvfCreate('wave',wave,'zcoeffs',zCoefs,'name',sprintf('human-%d',pupilMM));
        wvfP = wvfSet(wvfP,'calc pupil size',pupilMM);
        wvfP = wvfComputePSF(wvfP);
        % [u,p,f] = wvfPlot(wvfP,'2d psf space','um',550);
        % set(gca,'xlim',[-20 20],'ylim',[-20 20]);
        
        oi = wvf2oi(wvfP);
        oi = oiSet(oi,'name',sprintf('Human WVF %.1f mm',pupilMM));
        
    case {'diffractionlimited'}
        % Default optics is f# = 4, diffraction limited
        oi = oiSet(oi,'optics',optics);
        
        % Set up the default glass diffuser with a 2 micron blur circle,
        % but skipped
        oi = oiSet(oi, 'diffuserMethod','skip');
        oi = oiSet(oi, 'diffuserBlur', 2*10^-6);
        oi = oiSet(oi, 'consistency', 1);
        
    case {'uniformd65'}
        % Uniform, D65 optical image.  No cos4th falloff, huge field of
        % view (120 deg). Used in lux-sec SNR testing and scripting
        oi = oiCreateUniformD65;
        
    case {'uniformee', 'uniformeespecify'}
        % Uniform, equal energy optical image. No cos4th falloff. Might be used in
        % lux-sec SNR testing or scripting.  Not really used now
        % (5.3.2005).
        wave = 400:10:700; sz = 32;
        if length(varargin) >= 1, sz = varargin{1}; end
        if length(varargin) >= 2, wave = varargin{2}; end
        oi = oiCreateUniformEE(sz,wave);
               
    otherwise
        error('Unknown oiType');
end

% We may store an sRGB image of the irradiance here
% oi.sRGB = [];

if addObject
    if length(vcGetObjects('OPTICALIMAGE')) < val, vcAddAndSelectObject('OPTICALIMAGE',oi);
    else vcReplaceAndSelectObject(oi,val); end
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
scene = sceneCreate('uniformd65');
scene = sceneSet(scene,'hfov',120);
vcAddAndSelectObject(scene);

oi = oiCreate('default',[],[],0);
optics = oiGet(oi,'optics');
optics = opticsSet(optics,'offaxismethod','skip');
optics = opticsSet(optics,'otfmethod','skip');
oi = oiSet(oi,'optics',optics);

oi = oiCompute(scene,oi);


end

%---------------------------------------------
function oi = oiCreateUniformEE(sz,wave)
%
%  Create a uniform, equal energy image with a very large field of view.
%  The optical image is created without any cos4th fall off so it can be
%  used for lux-sec SNR testing.
%

scene = sceneCreate('uniformEESpecify',sz,wave);
scene = sceneSet(scene,'hfov',120);
vcAddAndSelectObject(scene);

oi = oiCreate('default',[],[],0);
optics = oiGet(oi,'optics');
optics = opticsSet(optics,'offaxismethod','skip');
optics = opticsSet(optics,'opticsModel','skip');
optics = opticsSet(optics,'otfmethod','skip');
oi = oiSet(oi,'optics',optics);

oi = oiCompute(scene,oi);

end
