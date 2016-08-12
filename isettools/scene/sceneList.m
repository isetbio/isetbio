% Scene types with parameters
%
%   sceneCreate('empty',wave)
%
% Macbeth scenes
%    sceneCreate('macbeth d65',patchSize);
%    sceneCreate('macbeth d50',patchSize);
%    sceneCreate('macbeth illC',patchSize);
%    sceneCreate('macbeth fluorescent',patchSize);
%    sceneCreate('macbeth tungsten',patchSize);
%    sceneCreate('macbeth ee_ir',patchSize);
%
%    sceneCreate('L star',barWidth,nBars,deltaEStep)
%
% Reflectance chart
%    sceneCreate('reflectance chart',pSize,sSamples,sFiles);
%    
% Monochromatic test
%    scene = sceneCreate('monochrome');     % Empty, but monochromatic
%    scene = sceneCreate('uniform monochromatic',wave,imsize);
%    scene = sceneCreate('multispectral');  % Empty, multispectral
%    scene = sceneCreate('rgb')
%
% Noise testing patterns
%    sceneCreate('linear Intensity Ramp', imSize)
%    sceneCreate('uniform Equal Energy', imSize, wave)
%    sceneCreate('uniform Equal Photon', imSize, wave)
%    sceneCreate('uniform d65', imSize)
%    sceneCreate('uniform bb', imSize,colorTemp, wave)
%    sceneCreate('white noise',row,col])
%
% General patterns
%    sceneCreate('rings rays',radialF,imsize)            
%    sceneCreate {'harmonic', paramStruct)
%    sceneCreate{'sweep frequency',imSize,maxFreq)
%
%    sceneCreate('line d65',imSize)
%    sceneCreate('line ee', imSize, offset, wave)
%    sceneCreate('line ep',imSize,offset)
%    sceneCreate('bar',imSize,barWidth);
%
%    sceneCreate('point array',imageSize,pixelsBetweenPoints);
%    sceneCreate('grid lines',imageSize,pixelsBetweenLines);
%    sceneCreate('radial lines', imageSizem spectralType, nLines);
%    sceneCreate('slanted edge',imageSize,edgeSlope);
%    sceneCreate('checkerboard',pixelsPerCheck,numberOfChecks)
%
%    sceneCreate('frequency orientation', params)
%      parms:  angles,freq,blockSize,contrast
%
%    sceneCreate('vernier',type, params) % See sceneVernier
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
%    sceneCreate('moire orient',imSize,edgeSlope);
%
% Text
%    font = fontCreate; d = displayCreate;
%    scene = sceneCreate('letter', font, d);
% 
%%

