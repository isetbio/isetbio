function [scene,parms] = sceneCreate(sceneName,varargin)
% Create a scene structure.
%
%  [scene,parms] = sceneCreate(sceneName,varargin)
%
% A scene describes the photons emitted from each visible point in the
% scene. Generally, we model planar objects, such as a screen display.  The
% scene is located at some distance from the center of the optics, has a
% field of view, and a spectral radiance distribution.  There are routines
% to handle depth as well that are partly implemented and under
% development.  We plan to integrate this aspect of the modeling with PBRT.
%
% A variety of scene types can be created automatically.  The routines that
% create these scenes, including this one, serve as a template for creating
% others you may wish to design.
%
% Scenes are represented as photons with 32 bits of precision by default.
% The spectral representation is 400:10:700 by default.  Both of these can
% be changed.
%
% To resample with respect to wavelength use the function
% sceneInterpolateW.
%
% The create a scene with 16 bits of precision, call sceneCreate with the
% full list of opticnal arguments and then append two arguments as in
%
%  sceneCreate(<complete list of optional arguments>,'bit depth',16)
%
% MACBETH COLOR AND LUMINANCE CHART
%
%   The default, scene = sceneCreate, is a Macbeth color checker illuminated
%   by a D65 light source with a mean luminance of 100 cd/m2.  The scene is
%   described only a small number of spatial 64x96 (row,col).  This can be
%   changed using the patchSize argument (default - 16 pixels).  The
%   wavelength  400:10:700 samples, making it efficient to use for experiments.
%
%   Example:
%    scene = sceneCreate('macbeth',32);
%
%    patchSize = 8;
%    spectrum.wave = (380:4:1068)';
%    scene = sceneCreate('macbethEE_IR',patchSize,spectrum);
%
%      {'macbethd65'}  - Create a Macbeth D65 image.  Optional
%         parameter of patch size (default = 16 pixels).
%      {'macbethd50'}         - D50 illuminant
%      {'macbethillc'}        - Illuminant C
%      {'macbethfluorescent'} - Fluorescent illuminant
%      {'macbethtungsten'}    - Tungsten illuminant
%      {'macbethEE_IR'}       - Equal energy extends out to the IR
%
%   You can use sceneAdjustIlluminant() to change the scene SPD.  This
%   function runs on MCC and any other scene with an illuminant.
%
%   The size of the individual patches and the wavelength sampling are both
%   parameters. They can be set using the calling procedure
%
%         patchSizePixels = 16;
%         spectrum.wave = [380:5:720];
%         scene = sceneCreate('macbethTungsten',patchSizePixels,spectrum);
%
%     {L*-steps} - Vertical bars spaced in equal L* steps (dark -> light)
%     scene = sceneCreate('LSteps',barWidth=20,nBars=10,deltaE=10);
%
% REFLECTANCE SAMPLE CHART
%
%   {'reflectance chart'}       - Specify random reflectance samples from
%                                 database. There is always a gray strip at
%                                 the right.   Uses sceneReflectanceChart
%
%       pSize = 24;          % Patch size in pixels
%       sSamples = [64 64];  % Surface samples from the files
%       sFiles{1} = fullfile(isetRootPath,'data','surfaces','reflectances','MunsellSamples_Vhrel.mat');
%       sFiles{2} = fullfile(isetRootPath,'data','surfaces','reflectances','Food_Vhrel.mat');
%       sceneCreate('reflectance chart',pSize,sSamples,sFiles);
%
% NARROWBAND COLOR PATCHES
%    wave = [600, 610];  sz = 64;
%    scene = sceneCreate('uniform monochromatic',wave,sz);
%
% SPATIAL TEST PATTERNS:
%
%      {'rings rays'}            - Resolution pattern
%      {'harmonic'}              - Harmonics (can be sums of harmonics)
%      {'sweep frequency'}       - Increasing frequency to the right,
%               increasing contrast upward
%      {'line d65'}              - Line with D65 energy spectrum
%      {'line ee'}               - Line with equal energy spectrum
%      {'bar ee'}                - Vertical bar, equal energy
%      {'point array'}           - Point array
%      {'gridlines'}             - Grid lines
%      {'checkerboard'}          - Checkerboard with equal photon spectrum
%      {'frequency orientation'} - Demosaicking test pattern, equal photon spectrum
%      {'slanted edge'} - Used for ISO spatial resolution, equal photon spectrum
%      {'moire orient'} - Circular Moire pattern
%      {'zone plate'}   - Circular zone plot, equal photon spectrum
%      {'star pattern'} - Radial lines used to test printers and displays
%
%  Additional parameters are available for several of the patterns.  For
%  example, the harmonic call can set the frequency, contrast, phase,
%  angle, row and col size of the harmonic.  The frequency unit in this
%  case is cycles/image.  To obtain cycles per degree, divide by the field
%  of view.
%
%        parms.freq = 1; parms.contrast = 1; parms.ph = 0;
%        parms.ang= 0; parms.row = 128; parms.col = 128;
%        parms.GaborFlag=0;
%        [scene,parms] = sceneCreate('harmonic',parms);
%
%  See the script s_sceneHarmonics for more examples.  In this example, the
%  illuminant is set so that the mean of the harmonic has a 20%
%  reflectance, like a typical gray card.
%
%  Many of the patterns can have an arbitrary image (row,col) size.  This
%  is possible for whitenoise, impulse1dee,lined65,
%
%         imSize = 128; lineOffset = 25;           % Plus is to the right
%         scene = sceneCreate('lined65',imSize);
%         scene = sceneCreate('line ee',imSize,lineOffset);
%         sceneCreate('bar',imageSize,width);
%
%  Other patterns have different parameters:
%
%         sceneCreate('slanted edge',imageSize,edgeSlope);
%         sceneCreate('checkerboard',pixelsPerCheck,numberOfChecks)
%         sceneCreate('grid lines',imageSize,pixelsBetweenLines);
%         sceneCreate('point array',imageSize,pixelsBetweenPoints);
%         sceneCreate('moire orient',imageSize,edgeSlope);
%         sceneCreate('vernier',imageSize,lineWidth,pixelOffset);
%
% TEXT
%      {'letter', 'font}- Scene created for certain character and display
%   For displays that have a psf, and a subset of all the possible letter
%   sizes and font combinations, we can create a scene.
%      letter = 'g'; fontSize = 18; 
%      display ='LCD-Apple'; fontName = 'Georgia';
%      scene = sceneCreate('letter', 'g', fontSize, fontName, display);
%
% NOISE ANALYSIS TEST PATTERNS
%
%      {'linearIntensityRamp'}  -
%      {'uniformEqualEnergy'}   - Equal energy
%      {'uniformEqualPhoton'}   - Equal photon density
%      {'uniform bb'}           - Blackbody, uniform
%      {'whitenoise'}           - Noise pattern for testing
%
%    The uniform patterns are small by default (32,32).  If you would like
%    them at a higher density (not much point), you can use
%
%        sceneCreate('uniform D65',256)
%        sceneCreate('uniform bb',128,6500)    - 6500 deg
%
%    where 256 is the image size in pixels.
%
% SCENES FROM IMAGE DATA
%   We also create scenes using data in image files.  It is also possible
%   to simply read a tiff or jpeg file and create a scene structure.  These
%   image-based scenes created by sceneFromFile.  See the comments there
%   for more information.
%
% See also:  sceneFromFile
%
% Copyright ImagEval Consultants, LLC, 2003.

if notDefined('sceneName'), sceneName = 'default'; end
parms = [];  % Returned in some cases, not many.

% Identify the object type
scene.type = 'scene';
scene = sceneSet(scene,'bit depth',32);   % Single precision
if length(varargin) > 1 && ischar(varargin{end-1})
    str = ieParamFormat(varargin{end-1});
    if isequal(str,'bitdepth')
        scene = sceneSet(scene,'bit depth',varargin{end});
    end
end
sceneName = ieParamFormat(sceneName);

switch sceneName
    case 'default'
        % The user can make a Macbeth with different patch sizes and
        % wavelength sampling, by calling with additional arguments, such
        % as scene = sceneCreate('macbethd65',16,spectrum);
        scene = sceneDefault(scene,'d65');
    case {'macbeth','macbethd65'}
        % sceneCreate('macbethD65',24);
        scene = sceneDefault(scene,'d65',varargin);
    case {'macbethd50'}
        scene = sceneDefault(scene,'d50',varargin);
    case {'macbethc','macbethillc'}
        scene = sceneDefault(scene,'c',varargin);
    case {'macbethfluorescent','macbethfluor'}
        scene = sceneDefault(scene,'fluorescent',varargin);
    case {'macbethtungsten','macbethtung'}
        scene = sceneDefault(scene,'tungsten',varargin);
    case {'macbethee_ir','macbethequalenergyinfrared'}
        % Equal energy illumination into the IR
        % The way to call this would be
        % patchSize = 16;
        % spectrum.wave = 380:4:1068;
        % scene = sceneCreate('macbethEE_IR',patchSize,spectrum)
        scene = sceneDefault(scene,'ir',varargin);
    case {'reflectancechart'}
        % sceneCreate('reflectance chart',pSize,sSamples,sFiles);
        % There is always a gray strip at the right.
        
        % Patch size in pixels
        pSize = 24;
        % Default surface files
        sFiles{1} = fullfile(isetRootPath,'data','surfaces', ...
                            'reflectances','MunsellSamples_Vhrel.mat');
        sFiles{2} = fullfile(isetRootPath,'data','surfaces', ...
                            'reflectances','Food_Vhrel.mat');
        % Surface samples from the files
        sSamples = [64 64];
        
        if isempty(varargin)
        else
            pSize = varargin{1};
            if length(varargin) > 1, sSamples = varargin{2}; end
            if length(varargin) > 2, sFiles = varargin{3}; end
        end
        scene = sceneReflectanceChart(sFiles,sSamples,pSize);
        
    case {'lstar'}
        scene = sceneSet(scene,'name','L-star');
        bWidth = 20; nBars = 10; deltaE = 10;
        if ~isempty(varargin)
            bWidth = varargin{1};
            if length(varargin)>1, nBars  = varargin{2}; end
            if length(varargin)>2, deltaE = varargin{3}; end
        end
        scene = sceneLstarSteps(scene,bWidth,nBars,deltaE);
        
        % Monochrome,RGB and multispectral add only a little.  Mostly created in sceneFromFile
    case {'monochrome','unispectral'}
        % Used for images with only one spectral band.
        scene = sceneSet(scene,'name','monochrome');
        scene = initDefaultSpectrum(scene,'monochrome');
    case {'multispectral','hyperspectral'}
        scene = sceneMultispectral(scene);
    case 'rgb'
        if isempty(varargin), scene = sceneRGB(scene);
        else scene = sceneRGB(varargin{1});end
        
    case {'mackay','rayimage','ringsrays'}
        % Also called the Siemens star pattern
        % radF = 24; imSize = 512;
        % vcAddAndSelectObject(sceneCreate('mackay',radF,imSize));
        % sceneWindow();
        radFreq = 8; sz = 256;
        if length(varargin) >= 1, radFreq = varargin{1}; end
        if length(varargin) >= 2, sz = varargin{2}; end
        if length(varargin) >= 3
            wave  = varargin{3};
            scene = sceneSet(scene,'wave',wave);
        end
        scene = sceneMackay(scene,radFreq,sz);
    case {'harmonic','sinusoid'}
        if isempty(varargin),
            [scene,parms] = sceneHarmonic(scene);
        elseif length(varargin) == 1
            parms = varargin{1};
            [scene,parms] = sceneHarmonic(scene,parms);
        else
            parms = varargin{1};
            wave = varargin{2};
            [scene,parms] = sceneHarmonic(scene,parms, wave);
        end
    case {'sweep','sweepfrequency'}
        % These are always equal photon type.  Could add a third argument
        % for spectral type.
        % sz = 512; maxF = sz/16; sceneCreate('sweepFrequency',sz,maxF);
        sz = 128; maxFreq = sz/16;
        if length(varargin) >= 1, sz = varargin{1}; end
        if length(varargin) >= 2, maxFreq = varargin{2}; end
        scene = sceneSweep(scene,sz,maxFreq);
    case {'ramp','linearintensityramp','rampequalphoton'}
        if isempty(varargin),  sz = 32;
        else                   sz = varargin{1};
        end
        scene = sceneRamp(scene,sz);
    case {'uniform','uniformee','uniformequalenergy'}   %Equal energy
        if isempty(varargin),  sz = 32;
        else                   sz = varargin{1};
        end
        scene = sceneUniform(scene,'equalenergy',sz);
        
    case {'uniformeespecify'}   % Equal energy, specify waveband
        % scene = sceneCreate('uniformEESpecify',sz,wavelength);
        sz = 32; wavelength = 400:10:700;
        if ~isempty(varargin), sz = varargin{1}; end
        if length(varargin) > 1, wavelength = varargin{2}; end
        scene = sceneSet(scene,'wave',wavelength(:));
        scene = sceneUniform(scene,'equalenergy',sz);
    case {'uniformequalphoton','uniformephoton'} %Equal photon density
        % sceneCreate('uniformEqualPhoton',128);
        if isempty(varargin)
            sz = 32;
            scene = sceneUniform(scene,'equal photons',sz);
        elseif length(varargin) == 1
            sz = varargin{1};
            scene = sceneUniform(scene,'equal photons',sz);
        else
            sz = varargin{1};
            wave = varargin{2};
            scene = sceneUniform(scene,'equal photons',sz, wave);
        end
    case 'uniformd65'
        % sceneCreate('uniformEqualPhoton',64);
        % We should include an option for wavelength so that we extend into
        % the IR
        if isempty(varargin),  sz = 32;
        else                   sz = varargin{1};
        end
        scene = sceneUniform(scene,'D65',sz);
    case {'uniformbb'}
        % scene = sceneCreate('uniform bb',64,5000,400:700);
        sz = 32; cTemp = 5000; wave = 400:10:700;
        if ~isempty(varargin),    sz    = varargin{1}; end
        if length(varargin) >= 2, cTemp = varargin{2}; end
        if length(varargin) >= 3, wave  = varargin{3}; end
        scene = sceneSet(scene,'wave',wave);
        scene = sceneUniform(scene,'BB',sz,cTemp);
    case {'uniformmonochromatic'}
        % scene = sceneCreate('uniform monochromatic',sz,wavelength);
        
        % Create a uniform, monochromatic image.  Used for color-matching
        % analyses.  Set the peak radiance in photons.
        sz = 128; wavelength = 500;
        if length(varargin) >= 1, wavelength = varargin{1}; end
        if length(varargin) >= 2, sz = varargin{2}; end
        
        scene = sceneSet(scene,'wave',wavelength);
        scene = sceneUniform(scene,'equalenergy',sz);
        scene = sceneSet(scene,'name','narrow band');
        
    case {'lined65','impulse1dd65'}
        if isempty(varargin), sz = 64;
        else sz = varargin{1};
        end
        scene = sceneLine(scene,'D65',sz);
    case {'lineee','impulse1dee'}
        % scene = sceneCreate('line ee',size,offset,wave);
        % size:   Image row/col
        % offset: Pixel offset from center (c + offset)
        % wave:   Wavelength samples
        % scene = sceneCreate('lineee',128,2);
        % scene = sceneCreate('lineee',128,2,380:4:1068);
        sz = 64; offset = 0;
        if length(varargin) >= 1, sz = varargin{1};     end
        if length(varargin) >= 2, offset = varargin{2}; end
        if length(varargin) == 3
            scene = sceneSet(scene,'wave',varargin{3});
        end
        scene = sceneLine(scene,'equalEnergy',sz,offset);
    case {'lineequalphoton','lineep'}
        % sceneCreate('line ep',sz,offset);
        sz = 64; offset = 0;
        if length(varargin) >= 1, sz = varargin{1};     end
        if length(varargin) >= 2, offset = varargin{2}; end
        scene = sceneLine(scene,'equalPhoton',sz,offset);
    case {'bar'}
        % sceneCreate('bar',sz,width)
        sz = 64; width = 3;
        if length(varargin) >=1, sz    = varargin{1};   end
        if length(varargin) >=2, width = varargin{2};   end
        scene = sceneBar(scene,sz,width);
    case {'vernier'}
        % sceneCreate('vernier', type, params)
        if ~isempty(varargin), type = varargin{1}; else type = object; end
        if length(varargin) > 1, 
            params = varargin{2};
        else
            params = [];
        end
        scene = sceneVernier(scene, type, params);
    case {'whitenoise','noise'}
        % sceneCreate('noise',[128 128])
        sz = 128; contrast = 20;
        if length(varargin) >= 1, sz = varargin{1}; end
        if length(varargin) >= 2, contrast = varargin{2}; end
        
        scene = sceneNoise(scene,sz,contrast);
        scene = sceneSet(scene,'name','white noise');
        
    case {'pointarray','manypoints'}
        % sceneCreate('pointArray',sz,spacing,spectralType);
        sz = 128; spacing = 16; spectralType = 'ep';
        if length(varargin) >= 1, sz           = varargin{1}; end
        if length(varargin) >= 2, spacing      = varargin{2}; end
        if length(varargin) >= 3, spectralType = varargin{3}; end
        scene = scenePointArray(scene,sz,spacing,spectralType);
        
    case {'gridlines','distortiongrid'}
        % sceneCreate('gridlines',sz,spacing,spectralType);
        sz = 128; spacing = 16; spectralType = 'ep';
        if length(varargin) >= 1, sz           = varargin{1}; end
        if length(varargin) >= 2, spacing      = varargin{2}; end
        if length(varargin) >= 3, spectralType = varargin{3}; end
        scene = sceneGridLines(scene,sz,spacing,spectralType);
        
    case {'checkerboard'}
        period = 16; spacing = 8; spectralType = 'ep';
        if length(varargin) >= 1, period       = varargin{1}; end
        if length(varargin) >= 2, spacing      = varargin{2}; end
        if length(varargin) >= 3, spectralType = varargin{3}; end
        scene = sceneCheckerboard(scene,period,spacing,spectralType);
        
    case {'demosaictarget','freqorientpattern', ...
            'frequencyorientation','freqorient'}
        %   parms.angles = linspace(0,pi/2,5);
        %   parms.freqs =  [1,2,4,8,16];
        %   parms.blockSize = 64;
        %   parms.contrast = .8;
        %   scene = sceneCreate('freqorient',parms);
        if isempty(varargin), scene = sceneFOTarget(scene);
        else
            % First argument is parms structure
            scene = sceneFOTarget(scene,varargin{1});
        end
        
    case {'moireorient'}
        %% Moire pattern test
        %   parms.angles = linspace(0,pi/2,5);
        %   parms.freqs =  [1,2,4,8,16];
        %   parms.blockSize = 64;
        %   parms.contrast = .8;
        % scene = sceneCreate('moire orient',parms);
        if isempty(varargin), scene = sceneMOTarget(scene);
        else
            % First argument is parms structure
            scene = sceneMOTarget(scene,varargin{1});
        end
    case {'slantedbar','iso12233','slantededge'}
        % scene = sceneCreate('slantedEdge',sz, slope, fieldOfView, wave);
        % scene = sceneCreate('slantedEdge',128,1.33);  % size, slope
        % scene = sceneCreate('slantedEdge',128,1.33,[], 380:4:780);
        barSlope = []; fov = []; wave = []; imSize = [];
        if length(varargin) >= 1, imSize = varargin{1}; end
        if length(varargin) >= 2, barSlope = varargin{2};  end
        if length(varargin) >= 3, fov = varargin{3}; end
        if length(varargin) >= 4, wave = varargin{4}; end
        scene = sceneSlantedBar(scene,imSize,barSlope,fov,wave);
        
    case {'zoneplate'}
        scene = sceneZonePlate(scene,384);
    case {'starpattern','radiallines'}
        % Thin radial lines - Useful for testing oriented blur
        %
        % scene = sceneCreate('starPattern');
        % scene = sceneCreate('starPattern',384);
        imSize = 256; spectralType = 'ep'; nLines = 8;
        if length(varargin) >=1, imSize = varargin{1}; end
        if length(varargin) >=2, spectralType = varargin{2}; end
        if length(varargin) >=3, nLines = varargin{3}; end
        scene = sceneRadialLines(scene,imSize,spectralType,nLines);
        
    case {'letter', 'font'}
        % Create scene of single letter
        %
        % scene = sceneCreate('letter',font,display);
        
        % Defaults, both have 96 dpi
        font = fontCreate; display = 'LCD-Apple';  
        
        % Assign arguments
        if ~isempty(varargin), font = varargin{1}; end
        if length(varargin) > 1, display = varargin{2}; end      
        if ischar(display), display = displayCreate(display); end
        
        scene = sceneFromFont(font,display);
        return;  % Do not adjust luminance or other properties
        
    otherwise
        error('Unknown scene format.');
end

%% Optional initializations
% The initializations below here are not required, but they are used by
% some of the default patterns above.  If you create a new scene type and
% don't want any of these initializations, call a return at the end of the
% case statement above.

% Initialize scene geometry, spatial sampling if not already set above
scene = sceneInitGeometry(scene);
scene = sceneInitSpatial(scene);

% Scenes are initialized to a mean luminance of 100 cd/m2.  The illuminant
% is adjusted so that dividing the radiance (in photons) by the illuminant
% (in photons) produces the appropriate peak reflectance (default = 1).
%
% Also, a best guess is made about one known reflectance.
if ~notDefined('scene.data.photons')
    if isempty(sceneGet(scene,'knownReflectance')) && ...
            checkfields(scene,'data','photons')       
        % If there is no illuminant yet, set the illuminant to equal
        % photons at 100 cd/m2
        wave = sceneGet(scene, 'wave');
        if isempty(sceneGet(scene,'illuminant'))
            il = illuminantCreate('equal photons', wave, 100);
            scene = sceneSet(scene, 'illuminant', il);
        end
        
        % There is no knownReflectance, so we set the peak radiance to a
        % reflectance of 0.9.
        v = sceneGet(scene, 'peakRadianceAndWave');
        idxWave = find(wave == v(2));
        p = sceneGet(scene, 'photons', v(2));
        [~, ij] = max2(p);
        v = [0.9 ij(1) ij(2) idxWave];
        scene = sceneSet(scene,'knownReflectance',v);
    end
    
    luminance = sceneCalculateLuminance(scene);
    scene = sceneSet(scene,'luminance',luminance);
    
    % This routine also adjusts the illumination level to be consistent
    % with the reflectance and scene photons.
    scene = sceneAdjustLuminance(scene,100);
    
    if ieSessionGet('gpu compute')
        p = sceneGet(scene, 'photons');
        scene = sceneSet(scene, 'photons', gpuArray(p));
    end
end

end

%---------------------------------------------------
function scene = sceneNoise(scene,sz,contrast)
%% Make a spatial white noise stimulus
% contrast is the standard deviation of the N(0,contrast) noise.
% The noise is shifted to a mean of 0.5, and the level is clipped to a
% minimum of 0.

if notDefined('sz'), sz = [128,128]; end
if notDefined('contrast'), contrast = 0.20;
elseif contrast > 1, contrast = contrast/100;
end

scene = initDefaultSpectrum(scene,'hyperspectral');
wave  = sceneGet(scene,'wave');
nWave = sceneGet(scene,'nwave');

% This is an image with reasonable dynamic range (10x).
d   = randn(sz)*contrast + 1; d = max(0,d);

% This is a D65 illuminant
il = illuminantCreate('d65',wave,100);
p  = illuminantGet(il,'photons');
scene = sceneSet(scene,'illuminant',il);

%
photons = zeros(sz(1),sz(2),nWave);
for ii=1:nWave, photons(:,:,ii) = d*p(ii); end

% Allocate space for the (compressed) photons
scene = sceneSet(scene,'photons',photons);

% By setting the fov here, we will not override the value in
% sceneInitSpatial() when this returns
scene = sceneSet(scene,'fov',1);

end


%----------------------------------
function scene = sceneDefault(scene,illuminantType,args)
%% Default scene is a Macbeth chart with D65 illuminant and patchSize 16
% pixels.

if notDefined('illuminantType'), illuminantType = 'd65'; end
if notDefined('args'), args = []; end

if (isempty(args) || isempty(args{1})), patchSize = 16;
else patchSize = args{1};
end

% Create the scene variable
scene = sceneSet(scene,'type','scene');
if isempty(args) || length(args) < 2 || isempty(args{2})
    scene = initDefaultSpectrum(scene,'hyperspectral');
else    scene = sceneSet(scene,'spectrum',args{2});
end
wave = sceneGet(scene,'wave');

switch lower(illuminantType)
    case 'd65'
        scene = sceneSet(scene,'name','Macbeth (D65)');
        lightSource = illuminantCreate('D65',wave,100);
    case 'd50'
        scene = sceneSet(scene,'name','Macbeth (D50)');
        lightSource = illuminantCreate('D50',wave,100);
    case 'fluorescent'
        scene = sceneSet(scene,'name','Macbeth (Fluorescent)');
        lightSource = illuminantCreate('Fluorescent',wave,100);
    case 'c'
        scene = sceneSet(scene,'name','Macbeth (Ill C)');
        lightSource = illuminantCreate('illuminantC',wave,100);
    case 'tungsten'
        scene = sceneSet(scene,'name','Macbeth (Tungsten)');
        lightSource = illuminantCreate('tungsten',wave,100);
    case 'ir'
        scene = sceneSet(scene,'name','Macbeth (IR)');
        lightSource = illuminantCreate('equalEnergy',wave,100);
    otherwise
        error('Unknown illuminant type.');
end

% Default distance in meters.
scene = sceneSet(scene,'distance',1.2);

% Scene magnification is always 1.
% Optical images have other magnifications that depend on the optics.
scene = sceneSet(scene,'magnification',1.0);

% The default patch size is 16x16.
spectrum = sceneGet(scene,'spectrum');
surface = macbethChartCreate(patchSize,(1:24),spectrum);

% scene = sceneCreateMacbeth(macbethChartObject,lightSource,scene);
iPhotons = illuminantGet(lightSource,'photons');
[surface,r,c] = RGB2XWFormat(surface.data);
sPhotons = surface*diag(iPhotons);
sPhotons = XW2RGBFormat(sPhotons,r,c);

% We compute the product of the surface reflectance and illuminant photons
% here
% scene = sceneSet(scene,'cphotons',surface.data .* photons);
scene = sceneSet(scene,'photons',sPhotons);

% Store the light source
scene = sceneSet(scene,'illuminant',lightSource);

end

%--------------------------------------------------
function scene = sceneMultispectral(scene)
%% Default multispectral structure

scene = sceneSet(scene,'name','multispectral');
scene = initDefaultSpectrum(scene,'multispectral');

end

%--------------------------------------------------
function scene = sceneRGB(scene)
%% Prepare a scene for RGB data.

if notDefined('scene'), scene.type = 'scene'; end

scene = sceneSet(scene,'name','rgb');
scene = sceneSet(scene,'type','scene');
scene = initDefaultSpectrum(scene,'hyperspectral');

% Set up an illuminant - but it is not nicely scaled.  And we don't have a
% known reflectance.
wave = sceneGet(scene,'wave');
il = illuminantCreate('d65',wave,100);
scene = sceneSet(scene,'illuminant',il);

end

%--------------------------------------------------
function scene = sceneMackay(scene,radFreq,sz)
%% Someone (I think Chris Tyler) told me the ring/ray pattern is also called
% the Mackay chart.
%   Reference from Joyce here:
%
% Some people call it the Siemens Star pattern (Wueller).
%
% We fill the central circle with a masking pattern.  The size of the
% central region is at the point when the rays would start to alias.  The
% the circumference of the central circle is 2*pi*r (with r in units of
% pixels).  When the radial frequency is f, we need a minimum of 2f pixels
% on the circumference.  So the circumference is 2*pi*r, so that we want
% the radius to be at least r = f/pi.  In practice that is too exact for
% the digital domain.  So we double the radius.
%

if notDefined('radFreq'), radFreq = 8; end
if notDefined('sz'),      sz = 256; end

scene = sceneSet(scene,'name','mackay');

if ~isfield(scene,'spectrum')
    scene = initDefaultSpectrum(scene,'hyperspectral');
end
nWave = sceneGet(scene,'nwave');

img = imgMackay(radFreq,sz);

% Insert central circle mask
r = round(2*radFreq/pi);  % Find the radius for the central circle

% Find the distance from the center of the image
[X,Y] = meshgrid(1:sz,1:sz); X = X - mean(X(:)); Y = Y - mean(Y(:));
d = sqrt(X.^2 + Y.^2);

% Everything with a distance less than 2r set to mean gray (128) for now.
l = (d < r);
img(l) = 128;  % figure; imagesc(img)

scene = sceneSet(scene,'cphotons',repmat(img,[1,1,nWave]));

% Set up an illuminant
wave = sceneGet(scene,'wave');
il = illuminantCreate('equal photons',wave,100);
scene = sceneSet(scene,'illuminant',il);

end

%--------------------------------------------------
function scene = sceneSweep(scene,sz,maxFreq)
%%  These are always equal photon

if notDefined('sz'), sz = 128; end
if notDefined('maxFreq'), maxFreq = sz/16; end

scene = sceneSet(scene,'name','sweep');
scene = initDefaultSpectrum(scene,'hyperspectral');
nWave = sceneGet(scene,'nwave');

img = imgSweep(sz,maxFreq);
img = img/max(img(:));

wave  = sceneGet(scene,'wave');
il    = illuminantCreate('equal photons',wave,100);
scene = sceneSet(scene,'illuminant',il);

img       = repmat(img,[1,1,nWave]);
[img,r,c] = RGB2XWFormat(img);
illP      = illuminantGet(il,'photons');
img       = img*diag(illP);
img       = XW2RGBFormat(img,r,c);
scene     = sceneSet(scene,'photons',img);

end

%--------------------------------------------------
function [scene,p] = sceneHarmonic(scene,params, wave)
%% Create a scene of a (windowed) harmonic function.
%
% Harmonic parameters are: parms.freq, parms.row, parms.col, parms.ang
% parms.ph, parms.contrast
%
% Missing default parameters are supplied by imageHarmonic.
%
% The frequency is with respect to the image (cyces/image).  To determine
% cycles/deg, use cpd: freq/sceneGet(scene,'fov');
%

scene = sceneSet(scene,'name','harmonic');

if notDefined('wave')
    scene = initDefaultSpectrum(scene,'hyperspectral');
else
    scene = initDefaultSpectrum(scene, 'custom',wave);
end

nWave = sceneGet(scene,'nwave');

% TODO: Adjust pass the parameters back from the imgHarmonic window. In
% other cases, they are simply attached to the global parameters in
% vcSESSION.  We can get them by a getappdata call in here, but not if we
% close the window as part of imageSetHarmonic
if notDefined('params')
    [h, params] = imageSetHarmonic; waitfor(h);
    img = imageHarmonic(params);
    p   = params;
else
    [img,p] = imageHarmonic(params);
end

% To reduce rounding error problems for large dynamic range, we set the
% lowest value to something slightly more than zero.  This is due to the
% ieCompressData scheme.
img(img==0) = 1e-4;
img   = img/(2*max(img(:)));    % Forces mean reflectance to 25% gray

% Mean illuminant at 100 cd
wave = sceneGet(scene,'wave');
il = illuminantCreate('equal photons',wave,100);
scene = sceneSet(scene,'illuminant',il);

img = repmat(img,[1,1,nWave]);
[img,r,c] = RGB2XWFormat(img);
illP = illuminantGet(il,'photons');
img = img*diag(illP);
img = XW2RGBFormat(img,r,c);
scene = sceneSet(scene,'photons',img);

% set scene field of view
scene = sceneSet(scene, 'h fov', 1);

end

%--------------------------------------------------
function scene = sceneRamp(scene,sz)
%% Intensity ramp (see L-star chart for L* steps)

if notDefined('sz'), sz = 128; end

scene = sceneSet(scene,'name','ramp');
scene = initDefaultSpectrum(scene,'hyperspectral');
nWave = sceneGet(scene,'nwave');
wave = sceneGet(scene,'wave');

img = imgRamp(sz);
img = img/(max(img(:)));

il = illuminantCreate('equal photons',wave,100);
scene = sceneSet(scene,'illuminant',il);

img = repmat(img,[1,1,nWave]);
[img,r,c] = RGB2XWFormat(img);
illP = illuminantGet(il,'photons');
img = img*diag(illP);
img = XW2RGBFormat(img,r,c);
scene = sceneSet(scene,'photons',img);

end

%--------------------------------------------------
function scene = sceneUniform(scene,spectralType,sz,varargin)
%% Create a spatially uniform scene.
%
% Various spd types are supported, including d65, blackbody, equal energy,
% equal photon
%

if notDefined('scene'), error('Scene required.'); end
if notDefined('spectralType'), spectralType = 'ep'; end
if notDefined('sz'), sz = [32 32]; end
if isscalar(sz), sz = [sz sz]; end
assert(numel(sz) == 2, 'sz should be 2 element vector');
sz = sz(:)'; % make sure sz is a row vector

scene = sceneSet(scene,'name',sprintf('uniform-%s',spectralType));

if isempty(varargin)
    if ~isfield(scene,'spectrum')
        scene = initDefaultSpectrum(scene,'hyperspectral');
    end
end
wave  = sceneGet(scene,'wave');
nWave = sceneGet(scene,'nwave');

% 100% reflectance
d = ones([sz nWave]);

spectralType = ieParamFormat(spectralType);
switch lower(spectralType)
    case {'d65','equalenergy','equalphotons','ee'}
        il = illuminantCreate(spectralType,wave);
    case {'blackbody','bb'}
        if isempty(varargin), cTemp = 5000;
        else                  cTemp = varargin{1};
        end
        il = illuminantCreate('blackbody',wave,cTemp);
    otherwise
        error('Unknown spectral type:%s\n',spectralType);
end

% Set illuminant
scene = sceneSet(scene,'illuminant',il);

% Create scene photons
illP = sceneGet(scene,'illuminant photons');
for ii=1:nWave, d(:,:,ii) = d(:,:,ii)*illP(ii); end

scene = sceneSet(scene,'cphotons',d);

end

%--------------------------------------------------
function scene = sceneLine(scene,spectralType,sz,offset)
%% Create a single line scene.
% This is used for computing linespreads and OTFs.

if notDefined('spectralType'), spectralType = 'ep'; end
if notDefined('sz'),     sz = 64; end
if notDefined('offset'), offset = 0; end

scene = sceneSet(scene,'name',sprintf('line-%s',spectralType));

if ~isfield(scene,'spectrum')
    scene = initDefaultSpectrum(scene,'hyperspectral');
end
wave    = sceneGet(scene,'wave');
nWave   = sceneGet(scene,'nwave');

% Black is more than zero to prevent HDR problem with ieCompressData
linePos = round(sz/2) + offset;
photons = ones(sz,sz,nWave)*1e-4;
photons(:,linePos,:) = 1;

spectralType = ieParamFormat(spectralType);
% Figure out a way to do this using sceneSet.
switch lower(spectralType)
    case {'ep','equalphotons','ephoton','equalphoton'}
        % Equal number of photons at every wavelength
        il = illuminantCreate('equal photons',wave);
    case {'ee','equalenergy','eenergy'}
        % Equal energy at every wavelength.  The large scale factor applied
        % to the number of photons is just to produce a reasonable energy
        % level.
        il = illuminantCreate('equal energy',wave);
        
    case 'd65'
        % D65 spectra for the line
        il = illuminantCreate('d65',wave);
        
    otherwise
        error('Unknown uniform field type %s.',spectralType);
end

scene = sceneSet(scene,'illuminant',il);
p     = sceneGet(scene,'illuminant photons');
for ii=1:nWave, photons(:,:,ii) = photons(:,:,ii)*p(ii); end

scene = sceneSet(scene,'photons',photons);

end

%--------------------------------------------------
function scene = sceneBar(scene,sz,width)
%% Create a single bar scene.
% This is used for computing the effect of scene dot density, say for a
% display with varying dots per inch.

if notDefined('sz'),     sz = 64; end
if notDefined('width'), width = 5; end

scene = sceneSet(scene,'name',sprintf('bar-%d',width));

if ~isfield(scene,'spectrum')
    scene = initDefaultSpectrum(scene,'hyperspectral');
end
wave    = sceneGet(scene,'wave');
nWave   = sceneGet(scene,'nwave');

% Black is more than zero to prevent HDR problem with ieCompressData
barPos = (1:width) + round((sz - width)/2);
photons = ones(sz,sz,nWave)*1e-8;   % Very dark reflectance mostly
photons(:,barPos,:) = 1;            % White reflectance in bar region

il = illuminantCreate('equal photons',wave);
scene = sceneSet(scene,'illuminant',il);
p     = sceneGet(scene,'illuminant photons');

% Create the radiance that matches reflectance and illuminant
for ii=1:nWave, photons(:,:,ii) = photons(:,:,ii)*p(ii); end

% Attach the photons to the scene
scene = sceneSet(scene,'photons',photons);

end

%------------------------
function scene = sceneVernier(scene, type, params)
%% scene for vernier acuity
%    type indicates what scene is created from, can take value from
%      'display' - scene is created from an image on display
%      'object'  - scene is creaetd from object with certain illuminance
%
%    params structure could include:
%      sceneSz    - scene resolution, default is 64
%      barWidth   - bar width in pixels
%      offset     - displacement in pixels
%      lineSpace  - spacing between the lines in pixels if type = 'display'
%      display    - display name or structure, useful if type = 'display'
%      barColor   - bar color, 0~1 RGB value for type = 'display'
%      bgColor    - background color, 0~1 RGB for type = 'display'
%      meanLum    - mean luminance
%      il         - illuminanece, for type = 'object'
%      barReflect - bar reflectance, for type = 'object'
%      bgReflect  - background reflectance, for type = 'object'
%
%

% check inputs
if notDefined('scene'), error('scene requried'); end
if notDefined('type'),  type = 'object'; end

% init parameters from params
if isfield(params, 'sceneSz'), sz = params.sceneSz; else sz = 64; end
if isfield(params, 'barWidth'), width = params.barWidth; else width = 0; end
if isfield(params, 'offset'), offset = params.offset; else offset = 1; end
if isfield(params, 'lineSpace'), lineSpace = params.lineSpace;
else lineSpace = inf; end

% Set scene parameters based on type
switch type
    case 'object'
        % Init bar and background reflectance parameter
        if isfield(params,'barReflect')
            lineReflectance = params.barReflect;
        else
            lineReflectance = 0.6;
        end
        if isfield(params, 'bgReflect')
            backReflectance = params.bgReflect;
        else
            backReflectance = 0.3;
        end
        
        scene = sceneSet(scene,'name',sprintf('vernier-%d',offset));
        
        %% We make the image square
        if isscalar(sz)
            r = sz; c = sz;
        else
            r = sz(1);  c = sz(2);
        end
        
        % Make the column number odd so we can really center the top line
        if ~isodd(c), c = c+1; end
        
        % Vernier line size and offset
        % Top and bottom half rows and columns
        % Columns containing top line, shifted offset/2
        topCols = (1:width) + round((c - width)/2) - floor(offset/2);
        
        % Columns containing bottom line, shifted offset from top columns
        % With this algorithm, the width of the
        botCols = topCols + offset;
        
        % Split the rows, too
        topHalf = round(r/2);
        topRows = 1:topHalf; botRows = (topHalf+1):r;
        
        %% Init spectrum
        
        if ~isfield(scene,'spectrum')
            scene = initDefaultSpectrum(scene,'hyperspectral');
        end
        wave    = sceneGet(scene,'wave');
        nWave   = sceneGet(scene,'nwave');
        
        %% Make the photon data
        if isfield(params,'il')
            il = params.il;
        else
            il    = illuminantCreate('equal photons',wave);
        end
        scene = sceneSet(scene,'illuminant',il);
        illP  = sceneGet(scene,'illuminant photons');
        
        photons = ones(r, c, nWave);
        for ii = 1 : nWave
            topBar = lineReflectance * illP(ii) * ...
                        photons(topRows, topCols, ii);
            botBar = lineReflectance * illP(ii) * ...
                        photons(botRows, botCols, ii);
            photons(:,:,ii) = backReflectance * photons(:,:,ii) * illP(ii);
            photons(topRows, topCols, ii)  = topBar;
            photons(botRows, botCols, ii)  = botBar;
        end
        
        scene = sceneSet(scene,'photons',photons);
    case 'display'
        % Init related parameters
        if isfield(params, 'display')
            display = params.display;
        else
            display = displayCreate('LCD-Apple');
        end
        if ischar(display), display = displayCreate(display); end
        
        if isfield(params, 'barColor')
            barColor = params.barColor;
        else
            barColor = 0.99;
        end
        if isscalar(barColor), barColor = repmat(barColor, [1 3]); end
        if isfield(params, 'bgColor')
            bgColor = params.bgColor;
        else
            bgColor = 0;
        end
        if isscalar(bgColor), bgColor = repmat(bgColor, [1 3]); end
        
        % Create image to be shown on display
        if isscalar(sz), sz = [sz sz]; end
        Img = repmat(reshape(bgColor,[1 1 3]),[sz 1]);
        cc = [round(sz(2)/2):-lineSpace:1 round(sz(2)/2):lineSpace:sz];
        width = width - 1;
        for jj = 1 : length(cc)
            barCols = max(round(cc(jj)-width/2),1) : ...
                      min(round(cc(jj)+width/2),sz(2));
            for ii = 1 : 3
                Img(:, barCols, ii) = barColor(ii);
            end
        end
        % Shift for offset
        Img(1:round(end/2),:,:) = circshift(Img(1:round(end/2),:,:), ...
            [0 offset 0]);
        
        % Create scene
        scene = sceneFromFile(Img, 'rgb',[], display);
    otherwise
        error('unknown vernier scene type');
end

if isfield(params, 'meanLum')
    scene = sceneAdjustLuminance(scene, params.meanLum);
end

end

%------------------------
function scene = sceneRadialLines(scene,imSize,spectralType,nLines)
%% sceneCreate('star pattern')
%  Create a Siemens Star (radial line) scene.
%
%   scene = sceneRadialLines(scene,imSize,spectralType,nLines)
%
% In this test chart the intensities along lines from the center are
% constant. Measuring on a circle around the center the intensity is a
% harmonic. Hence, frequency varies as a function of radial distance.
%
% Reference:
%   Dieter Wueller thinks this pattern is cool.
%   Digital camera resolution measurement using sinusoidal Siemens stars
%   Proc. SPIE, Vol. 6502, 65020N (2007); doi:10.1117/12.703817
%
% Examples:
%  scene = sceneCreate('radialLines');
%

if notDefined('scene'), error('Scene must be defined'); end
if notDefined('spectralType'), spectralType = 'ep'; end
if notDefined('imSize'), imSize = 256; end
if notDefined('nLines'), nLines = 8; end

scene = sceneSet(scene,'name',sprintf('radialLine-%s',spectralType));
scene = initDefaultSpectrum(scene,'hyperspectral');

% Determine the line angles
radians = pi*(0:(nLines-1))/nLines;
endPoints = zeros(nLines,2);
for ii=1:nLines
    endPoints(ii,:) = round([cos(radians(ii)),sin(radians(ii))]*imSize/2);
end
% plot(endPoints(:,1),endPoints(:,2),'o')

img = zeros(imSize,imSize);

% The routine for drawing lines could be better.
for ii=1:nLines
    x = endPoints(ii,1); y = endPoints(ii,2);
    u = -x; v = -y;
    % Flip so x is the lower one
    if x > 0,
        tmp = [x,y]; x = u; y = v; u = tmp(1); v = tmp(2);
    end
    
    if ~isequal(u,x), slope = (y - v) / (u - x);
        for jj=x:0.2:u,
            kk = round(jj*slope);
            img(round(kk + (imSize/2)) + 1, round(jj + (imSize/2)) + 1) = 1;
        end
    else img(:, (imSize/2) + 1) = 1;
    end
end

img = img(1:imSize,1:imSize);
% To reduce rounding error problems for large dynamic range, we set the
% lowest value to something slightly more than zero.  This is due to the
% ieCompressData scheme.
img(img==0) = 1e-4;
img = img/max(img(:));
% figure; imagesc(img)

% Create the photon image
wave    = sceneGet(scene,'wave');
nWave   = sceneGet(scene,'nwave');
photons = zeros(imSize,imSize,nWave);

% Figure out a way to do this using sceneSet.
spectralType = ieParamFormat(spectralType);
switch spectralType
    case {'ep','equalphoton','ephoton'}
        % Equal number of photons at every wavelength
        il = illuminantCreate('equal photons',wave);
    case {'ee','equalenergy','eenergy'}
        il = illuminantCreate('equal energy',wave);
    case 'd65'
        % D65 spectra for the line
        il = illuminantCreate('d65',wave);
    otherwise
        error('Unknown uniform field type %s.\n',spectralType);
end

scene = sceneSet(scene,'illuminant',il);
p     = sceneGet(scene,'illuminant photons');
for ii=1:nWave, photons(:,:,ii) = img(:,:)*p(ii); end

scene = sceneSet(scene,'photons',photons);

end

%-----------------------
function scene = sceneFOTarget(scene,parms)
%% Frequency/Orientation target

if notDefined('parms'), parms = []; end

scene = sceneSet(scene,'name','FOTarget');
scene = initDefaultSpectrum(scene,'hyperspectral');
nWave = sceneGet(scene,'nwave');

img = FOTarget('sine',parms);

% Prevent dynamic range problem with ieCompressData
img = ieClip(img,1e-4,1);
img = img/max(img(:));


% Create the illuminant
il = illuminantCreate('equal photons',sceneGet(scene,'wave'));
scene = sceneSet(scene,'illuminant',il);

% This routine returns an RGB image.  We base the final image on just the
% green channel
img = repmat(img(:,:,2),[1,1,nWave]);
[img,r,c] = RGB2XWFormat(img);
illP = illuminantGet(il,'photons');
img = img*diag(illP);
img = XW2RGBFormat(img,r,c);
scene = sceneSet(scene,'photons',img);

end
%-----------------------
function scene = sceneMOTarget(scene,parms)
%% Moire/Orientation target

if notDefined('parms'), parms = []; end

scene = sceneSet(scene,'name','MOTarget');
scene = initDefaultSpectrum(scene,'hyperspectral');
nWave = sceneGet(scene,'nwave');

% Select one among sinusoidalim, squareim, sinusoidalim_line,
% squareim_line, flat img = MOTarget('squareim',parms);
img = MOTarget('sinusoidalim',parms);


% Prevent dynamic range problem with ieCompressData
img = ieClip(img,1e-4,1);

% This routine returns an RGB image.  We take the green channel and expand
% it
scene = sceneSet(scene,'cphotons',repmat(img(:,:,2),[1,1,nWave]));

%
wave = sceneGet(scene,'wave');
illPhotons = ones(size(wave))*sceneGet(scene,'data max');
scene = sceneSet(scene,'illuminantPhotons',illPhotons);

end

%-------------------
function scene = sceneCheckerboard(scene,checkPeriod,nCheckPairs,spectralType)
%% Checkerboard

if notDefined('scene'), error('Scene required'); end
if notDefined('checkPeriod'), checkPeriod = 16; end
if notDefined('nCheckPairs'), nCheckPairs = 8; end
if notDefined('spectralType'), spectralType = 'ep'; end

scene = sceneSet(scene,'name',sprintf('Checker-%s',spectralType));
scene = initDefaultSpectrum(scene,'hyperspectral');
wave  = sceneGet(scene,'wave');
nWave = sceneGet(scene,'nwave');

% Changed to Matlab checkerboard in image processing toolbox.  Forced the
% output to double.
d = checkerboard(checkPeriod,nCheckPairs); d = double((d > 0.5));
d = d/max(d(:));

% Prevent ieCompressData problem.
d = ieClip(d,1e-6,1);
spectralType = ieParamFormat(spectralType);
switch spectralType
    case {'d65'}
        il = illuminantCreate('d65',wave);
    case {'ee','equalenergy'}
        il = illuminantCreate('equalenergy',wave);
    case {'ep','equalphoton','equalphotons'}
        il = illuminantCreate('equal photons',wave);
    otherwise
        error('Unknown spectral type:%s\n',spectralType);
end

img = zeros(size(d,1),size(d,2),nWave);
illP = illuminantGet(il,'photons');
for ii=1:nWave, img(:,:,ii) = d*illP(ii); end
scene = sceneSet(scene,'photons',img);

end

%---------------------------------------------------------------
function scene = sceneSlantedBar(scene,imSize,barSlope,fieldOfView,wave)
%%
%  Slanted bar, 2 deg field of view
%  Slope 2.6 (upper left to lower right)
%  Default size:  384
%
% The scene is set to equal photons across wavelength.

if notDefined('imSize'),      imSize = 384; end
if notDefined('barSlope'),    barSlope = 2.6; end
if notDefined('fieldOfView'), fieldOfView = 2; end
if notDefined('wave'),        wave = 400:10:700; end
scene = sceneSet(scene,'name','slantedBar');

scene = sceneSet(scene,'wave',wave);

wave = sceneGet(scene,'wave');
nWave  = sceneGet(scene,'nwave');

% Make the image
imSize = round(imSize/2);
[X,Y] = meshgrid(-imSize:imSize,-imSize:imSize);
img = zeros(size(X));

%  y = barSlope*x defines the line.  We find all the Y values that are
%  above the line
list = (Y > barSlope*X );

% We assume target is perfectly reflective (white), so the illuminant is
% the equal energy illuminant; that is, the SPD is all due to the
% illuminant
img( list ) = 1;

% Prevent dynamic range problem with ieCompressData
img = ieClip(img,1e-6,1);

% Now, create the illuminant
il = illuminantCreate('equal energy',wave);
scene = sceneSet(scene,'illuminant',il);
illP = illuminantGet(il,'photons');

% Create the scene photons
photons = zeros(size(img,1),size(img,2),nWave);
for ii=1:nWave, photons(:,:,ii) = img*illP(ii); end
scene = sceneSet(scene,'photons',photons);

% Set the field of view
scene = sceneSet(scene,'horizontalfieldofview',fieldOfView);

end

%-----------------------
function scene = sceneZonePlate(scene,imSize,fieldOfView)
%% Circular zone plate image
%

if notDefined('imSize'), imSize = 256; end
if notDefined('fieldOfView'), fieldOfView = 4; end

scene = sceneSet(scene,'name','zonePlate');
scene = initDefaultSpectrum(scene,'hyperspectral');
nWave = sceneGet(scene,'nwave');

img = imgZonePlate(imSize);
% Prevent dynamic range problem with ieCompressData
img = ieClip(img,1e-4,1);

scene = sceneSet(scene,'cphotons',repmat(img,[1,1,nWave]));
scene = sceneSet(scene,'horizontalfieldofview',fieldOfView);

end

%-----------------------
function scene = sceneLstarSteps(scene,barWidth,nBars,deltaE)
%% Scene with vertical bars in equal L* steps
%
% scene = sceneCreate('lstar',50,5,10);
% vcAddAndSelectObject(scene); sceneWindow;

scene = initDefaultSpectrum(scene,'hyperspectral');

% Create the Y values that will define the intensities of the spd.  First,
% equal spaced L* values
L = (0:(nBars-1))*deltaE;
LAB = zeros(nBars,3);
LAB(:,1) = L(:);

% Transform them to Y values
C = makecform('lab2xyz');
XYZ = applyCform(LAB,C);
Y = XYZ(:,2); Y = Y/max(Y(:));
% vcNewGraphWin; plot(Y)

% Create equal photons illuminant
il = illuminantCreate('equal photons',sceneGet(scene,'nwave'),100);
scene = sceneSet(scene,'illuminant',il);

% Now, make the photon image
photons = ones(128,barWidth*nBars,nWave);
for ii=1:nBars
    start = barWidth*(ii-1) + 1; stop = barWidth*ii;
    for jj=1:nWave
        photons(:,start:stop,jj) = Y(ii)*illPhotons(jj);
    end
end

% The level is scaled on return to a mean of 100 cd/m2.
scene = sceneSet(scene,'photons',photons);

end
