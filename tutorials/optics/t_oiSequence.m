% t_oiSequence
%
% Demonstrate usage of oiSequence
%
%
% NPC, ISETBIO TEAM, 2016
%

meanLuminance = 30;
uniformScene = sceneCreate('uniform equal photon', 128);
% square scene with desired FOV
FOV = 2.0;
uniformScene = sceneSet(uniformScene, 'wAngular', FOV);
% 1 meter away
uniformScene = sceneSet(uniformScene, 'distance', 1.0);
uniformScene = sceneAdjustLuminance(uniformScene, meanLuminance);
 
% generating stimulus modulation function
stimulusSamplingInterval = 60/1000;
oiTimeAxis = -0.6:stimulusSamplingInterval:0.6;
stimulusRampTau = 0.165;
modulationFunction = exp(-0.5*(oiTimeAxis/stimulusRampTau).^2);

% Default human optics
oi = oiCreate('human');
oi = oiCompute(oi, uniformScene);

% Compute the background and the modulated optical images
oiBackground = oiCompute(oi, uniformScene);
oiModulated  = oiBackground;
  
% Instantiate an oiSequence object for computing the sequence of the (full field) modulated ois
theOIsequence = oiSequence(oiBackground, oiModulated, modulationFunction);

% Instantiate another oiSequence object for computing the sequence of the (windowed, radius = 250 microns) modulated ois
modulationRegion.radiusInMicrons = 250;
theOIsequence2 = oiSequence(oiBackground, oiModulated, modulationFunction, 'modulationRegion', modulationRegion);

% Plot the modulation function
figure(1);
plot(1:theOIsequence.length, modulationFunction, 'ks-', 'LineWidth', 1.5);
set(gca, 'XLim', [1 theOIsequence.length]);
xlabel('frame index');
ylabel('modulation');

% Plot the 2 oisequences
figure(2); clf;
for oiIndex = 1:theOIsequence.length
    % Ask theOIsequence to return the oiIndex-th frame
    currentOI = theOIsequence.frameAtIndex(oiIndex);
    
    % plot it
    subplot(4,round(theOIsequence.length/2), oiIndex);
    rgbImage = xyz2rgb(oiGet(currentOI, 'xyz'));
    imagesc(rgbImage, [0 1]);
    title(sprintf('frame %d', oiIndex));
    axis 'image'
    set(gca, 'XTick', [], 'YTick', []);
    
    % Ask theOIsequence2 to return the oiIndex-th frame
    currentOI = theOIsequence2.frameAtIndex(oiIndex);
    
    % plot it
    subplot(4,round(theOIsequence.length/2), oiIndex+2*(round(theOIsequence.length/2)));
    rgbImage = xyz2rgb(oiGet(currentOI, 'xyz'));
    imagesc(rgbImage, [0 1]);
    title(sprintf('frame %d', oiIndex));
    axis 'image'
    set(gca, 'XTick', [], 'YTick', []);
    
end
