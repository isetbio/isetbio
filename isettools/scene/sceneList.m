%% Scene types with parameters
%
% General
%     wave = 400:10:700;
%     sceneCreate('empty', wave)
%
%%
% Macbeth scenes
%     patchSize = 8;
%     barWidth = 20;
%     nBars=10;
%     deltaE=10;
% 
%     sceneCreate('macbeth d65', patchSize);
%     sceneCreate('macbeth d50', patchSize);
%     sceneCreate('macbeth illC', patchSize);
%     sceneCreate('macbeth fluorescent', patchSize);
%     sceneCreate('macbeth tungsten', patchSize);
%     sceneCreate('macbeth ee_ir', patchSize);
% 
%     % This one isn't working.
%     sceneCreate('L star', barWidth, nBars, deltaE)
%     sceneShowImage(scene)
%
%%
% Reflectance chart
%     pSize = 24;
%     sSamples = [64 64];
%     sFiles{1} = fullfile(isetbioDataPath, 'surfaces', 'reflectances', ...
%         'MunsellSamples_Vhrel.mat');
%     sFiles{2} = fullfile(isetbioDataPath, 'surfaces', 'reflectances', ...
%         'Food_Vhrel.mat');
%     scene = sceneCreate('reflectance chart', pSize, sSamples, sFiles);
%     sceneShowImage(scene)
%
%%
% Monochromatic test
%     imSize = 128;
%     wave = 400:10:700;
%     scene = sceneCreate('monochrome');     % Empty, but monochromatic
%     scene = sceneCreate('uniform monochromatic', wave, imSize);
%     scene = sceneCreate('multispectral');  % Empty, multispectral
%     scene = sceneCreate('rgb')
%
%%
% Noise testing patterns
%     imSize = 128;
%     wave = 400:10:700;
%     colorTemp = 'hot'
% 
%     sceneCreate('linear Intensity Ramp', imSize)
%     sceneCreate('uniform Equal Energy', imSize, wave)
% 
%     % This one isn't working
%     sceneCreate('uniform Equal Photon', imSize, wave)
% 
%     sceneCreate('uniform d65', imSize)
% 
%     % This one isn't working
%     sceneCreate('uniform bb', imSize, colorTemp, wave)
% 
%     sceneCreate('white noise', 200, 300])
%
%%
% General patterns
%     sceneCreate('rings rays', radialF, imsize)            
%     sceneCreate {'harmonic', paramStruct)
%     sceneCreate{'sweep frequency', imSize, maxFreq)
% 
%     sceneCreate('line d65', imSize)
%     sceneCreate('line ee', imSize, offset, wave)
%     sceneCreate('line ep', imSize, offset)
%     sceneCreate('bar', imSize, barWidth);
% 
%     sceneCreate('point array', imageSize, pixelsBetweenPoints);
%     sceneCreate('grid lines', imageSize, pixelsBetweenLines);
%     sceneCreate('radial lines', imageSizem spectralType, nLines);
%     sceneCreate('slanted edge', imageSize, edgeSlope);
%     sceneCreate('checkerboard', pixelsPerCheck, numberOfChecks)
%
%    sceneCreate('frequency orientation', params)
%      parms:  angles, freq, blockSize, contrast
%
%    sceneCreate('vernier', type, params) % See sceneVernier
%     params:
%      sceneSz    - scene resolution, default is 64
%      barWidth   - bar width in pixels
%      offset     - displacement in pixels
%      meanLum    - mean luminance
%     if type  = 'display'
%      lineSpace  - spacing between the lines in pixels 
%      display    - display name or structure
%      barLength  - length of the line segment in number of pixels
%      barColor   - bar color, 0~1 RGB value 
%      bgColor    - background color, 0~1 RGB 
%     if type = 'object'
%      il         - illuminance
%      barReflect - bar reflectance
%      bgReflect  - background reflectance
%
%    sceneCreate('zone plate', imSize);
%    sceneCreate('moire orient', imSize, edgeSlope);
%
%% 
% Text
%    font = fontCreate;
%    d = displayCreate;
%    scene = sceneCreate('letter', font, d);
% 
%%
