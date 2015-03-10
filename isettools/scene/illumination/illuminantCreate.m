function il = illuminantCreate(ilName,wave, varargin)
% Create an illuminant (light source) structure.  
%
%  il = illuminantCreate(ilName,temperature,luminance,spectrum)
%
% Historically the illuminant is just a spectral function. Several standard
% light sources are supported at present.  These are d65, d50, tungsten,
% fluorescene, blackbody's (over a range of color temperatures), 550nm. See
% the internal routine, illuminantRead (below).
%
% We are now starting to experiment with spatial spectrum illuminants, that
% is a separate illuminant for each point in the scene data.
%
% The illuminant data are stored in units of [photons/(sr m^2 nm)]
%
% Illuminants
%    blackbody tempDegKelvin
%    tungsten
%    d65
%    fluorescent
%    555nm
%    equal energy
%    illuminant c
%
% Examples:
%   il = illuminantCreate('d65')
%   il = illuminantCreate('blackbody',400:10:700, 3500,100)
%   il = illuminantCreate('blackbody',,[],6500,100)
%   il = illuminantCreate('illuminant c',400:1:700,500)
%
%   spectrum.wave = (380:4:1068);
%   il = illuminantCreate('equalEnergy',[],100,spectrum)
%
% See also:  illuminantSet/Get, s_sceneIlluminant, s_sceneIlluminantSpace,
%            illuminantRead
%
% Copyright ImagEval Consultants, LLC, 2005.

%% Initialize parameters
if notDefined('ilName'), ilName = 'd65'; end

il.name = ilName;
il.type = 'illuminant';
il = initDefaultSpectrum(il,'hyperspectral');
if exist('wave','var') && ~isempty(wave), il.spectrum.wave = wave; end

%% There is no default
% The absence of a default could be a problem.

switch ieParamFormat(ilName)
    case {'d65','d50','tungsten','fluorescent', ...
          '555nm','equalenergy','illuminantc','equalphotons'}
        % illuminantCreate('d65',luminance)
        illP.name = ilName;
        illP.luminance = 100;
        illP.spectrum.wave = illuminantGet(il,'wave');
        if ~isempty(varargin), illP.luminance = varargin{1}; end;
        
        iEnergy = illuminantRead(illP);		    % [W/(sr m^2 nm)]
        iPhotons = Energy2Quanta(illuminantGet(il,'wave'),iEnergy); % Check this step
        il = illuminantSet(il,'name',illP.name);

    case 'blackbody'
        % illuminantCreate('blackbody',5000,luminance);
        illP.name = 'blackbody';
        illP.temperature = 5000;
        illP.luminance   = 100;
        illP.spectrum.wave = illuminantGet(il,'wave');
        
        if ~isempty(varargin),   illP.temperature = varargin{1};  end
        if length(varargin) > 1, illP.luminance = varargin{2}; end;
        
        iEnergy = illuminantRead(illP);		    % [W/(sr m^2 nm)]
        iPhotons = Energy2Quanta(illuminantGet(il,'wave'),iEnergy); % Check this step
        
        il = illuminantSet(il,'name',sprintf('blackbody-%.0f',illP.temperature));
        
    otherwise
        error('unknown illuminant type %s\n',ilName);
end

%% Set the photons and return
il = illuminantSet(il,'photons',iPhotons);  % [photons/(s sr m^2 nm)]

end