% t_oiSequence
%
% Demonstrate usage of oiSequence
%
%
% NPC, ISETBIO TEAM, 2016
%

% Generate a uniform scene
meanLuminance = 10;
uniformScene = sceneCreate('uniform equal photon', 128);
% square scene with desired FOV
FOV = 2.0;
uniformScene = sceneSet(uniformScene, 'wAngular', FOV);
% 1 meter away
uniformScene = sceneSet(uniformScene, 'distance', 1.0);
uniformScene = sceneAdjustLuminance(uniformScene, meanLuminance);

% Generate a gabor scene
gaborParams = struct(...
    'freq', 2, ...
    'contrast', 1.0, ...
    'ph', 0, ...
	'ang',  0, ...
    'row', 128, ...
    'col', 128, ...
	'GaborFlag', false);
gaborScene = sceneCreate('harmonic', gaborParams);
gaborScene = sceneSet(gaborScene, 'wAngular', FOV);
gaborScene = sceneSet(gaborScene, 'distance', 1.0);
gaborScene = sceneAdjustLuminance(gaborScene, meanLuminance);

% generating stimulus modulation functions
stimulusSamplingInterval = 60/1000;
oiTimeAxis = -0.6:stimulusSamplingInterval:0.6;
stimulusRampTau = 0.165;
% monophasic modulation function
modulationFunction1 = 0.7*exp(-0.5*(oiTimeAxis/stimulusRampTau).^2);
modulationFunction2 = -modulationFunction1;
% biphasic modulation function
modulationFunction3 = 0.7*(exp(-0.5*((oiTimeAxis-0.2)/stimulusRampTau).^2) - exp(-0.5*((oiTimeAxis+0.2)/stimulusRampTau).^2));

% Default human optics
oi = oiCreate('human');
oi = oiCompute(oi, uniformScene);

% Compute the background and the modulated optical images
oiBackground = oiCompute(oi, uniformScene);
oiModulated  = oiBackground;
oiModulatedGabor = oiCompute(oi, gaborScene);
modulationRegion.radiusInMicrons = 250;

% oiSequence object for computing a sequence of ois where the oiModulated
% (uniform field) is ADDED to the oiBackground over an 250 micron radius
% region using a monophasic modulation function
theOIsequence = oiSequence(oiBackground, oiModulated, modulationFunction1, 'composition', 'add', 'modulationRegion', modulationRegion);

% oiSequence object for computing a sequence of ois where the oiModulated
% (a grating) is BLENDED with the oiBackground using a biphasic modulation function
theOIsequence2 = oiSequence(oiBackground, oiModulated, modulationFunction2, 'composition', 'add', 'modulationRegion', modulationRegion);

% oiSequence object for computing a sequence of ois where the oiModulated
% (a grating) is ADDED with the oiBackground using a biphasic modulation function
theOIsequence3 = oiSequence(oiBackground, oiModulatedGabor, modulationFunction1,  'composition', 'add');

% oiSequence object for computing a sequence of ois where the oiModulated
% (a grating) is BLENDED with the oiBackground using a biphasic modulation function
theOIsequence4 = oiSequence(oiBackground, oiModulatedGabor, modulationFunction3, 'composition', 'blend');

% Plot the oisequences
hFig = figure(1); clf;
set(hFig, 'Color', [1 1 1], 'Position', [10 10 1500 950]);
for oiIndex = 1:theOIsequence.length
    
    % Plot the modulation function
    subplot(8,round(theOIsequence.length/2)+1, 1);
    plot(1:theOIsequence.length, theOIsequence.modulationFunction, 'ks-', 'LineWidth', 1.5);
    set(gca, 'XLim', [1 theOIsequence.length]);
    title(sprintf('mode: %s', theOIsequence.composition));
    xlabel('frame index');
    ylabel('modulation');
    
    % Ask theOIsequence to return the oiIndex-th frame
    currentOI = theOIsequence.frameAtIndex(oiIndex);
    
    % plot it
    subplot(8,round(theOIsequence.length/2)+1, 1+oiIndex);
    rgbImage = lrgb2srgb(xyz2rgb(oiGet(currentOI, 'xyz')));
    imagesc(rgbImage, [0 1]);
    title(sprintf('frame %d', oiIndex));
    axis 'image'
    set(gca, 'XTick', [], 'YTick', []);
    
    
    % Plot the modulation function
    subplot(8,round(theOIsequence2.length/2)+1, 1 + 2*(1+round(theOIsequence2.length/2)));
    plot(1:theOIsequence2.length, theOIsequence2.modulationFunction, 'ks-', 'LineWidth', 1.5);
    set(gca, 'XLim', [1 theOIsequence.length]);
    title(sprintf('mode: %s', theOIsequence2.composition));
    xlabel('frame index');
    ylabel('modulation');
    
    % Ask theOIsequence2 to return the oiIndex-th frame
    currentOI = theOIsequence2.frameAtIndex(oiIndex);
    
    % plot it
    subplot(8,round(theOIsequence.length/2)+1, 1+oiIndex+2*(1+round(theOIsequence.length/2)));
    rgbImage = lrgb2srgb(xyz2rgb(oiGet(currentOI, 'xyz')));
    imagesc(rgbImage, [0 1]);
    title(sprintf('frame %d', oiIndex));
    axis 'image'
    set(gca, 'XTick', [], 'YTick', []);
    
    
    % Plot the modulation function
    subplot(8,round(theOIsequence3.length/2)+1, 1 + 4*(1+round(theOIsequence3.length/2)));
    plot(1:theOIsequence3.length, theOIsequence3.modulationFunction, 'ks-', 'LineWidth', 1.5);
    set(gca, 'XLim', [1 theOIsequence3.length]);
    title(sprintf('mode: %s', theOIsequence3.composition));
    xlabel('frame index');
    ylabel('modulation');
    
    % Ask theOIsequence3 to return the oiIndex-th frame
    currentOI = theOIsequence3.frameAtIndex(oiIndex);
    
    % plot it
    subplot(8,round(theOIsequence3.length/2)+1, 1+oiIndex+4*(1+round(theOIsequence3.length/2)));
    rgbImage = lrgb2srgb(xyz2rgb(oiGet(currentOI, 'xyz')));
    imagesc(rgbImage, [0 1]);
    title(sprintf('frame %d', oiIndex));
    axis 'image'
    set(gca, 'XTick', [], 'YTick', []);
    
    
    % Plot the modulation function
    subplot(8,round(theOIsequence4.length/2)+1, 1 + 6*(1+round(theOIsequence4.length/2)));
    plot(1:theOIsequence.length, theOIsequence4.modulationFunction, 'ks-', 'LineWidth', 1.5);
    set(gca, 'XLim', [1 theOIsequence4.length]);
    title(sprintf('mode: %s', theOIsequence4.composition));
    xlabel('frame index');
    ylabel('modulation');
    
    % Ask theOIsequence4 to return the oiIndex-th frame
    currentOI = theOIsequence4.frameAtIndex(oiIndex);
    
    % plot it
    subplot(8,round(theOIsequence4.length/2)+1, 1+oiIndex+6*(1+round(theOIsequence4.length/2)));
    rgbImage = lrgb2srgb(xyz2rgb(oiGet(currentOI, 'xyz')));
    imagesc(rgbImage, [0 1]);
    title(sprintf('frame %d', oiIndex));
    axis 'image'
    set(gca, 'XTick', [], 'YTick', []);
    
    
end
