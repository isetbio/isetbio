% t_oiSequence
%
% Demonstrate usage of oiSequence
%
%
% NPC, ISETBIO TEAM, 2016
%

% Generate a uniform scene
meanLuminance = 40;
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
theOIsequence(1) = oiSequence(oiBackground, oiModulated, modulationFunction1, 'composition', 'add', 'modulationRegion', modulationRegion);

% oiSequence object for computing a sequence of ois where the oiModulated
% (a grating) is BLENDED with the oiBackground using a biphasic modulation function
theOIsequence(2) = oiSequence(oiBackground, oiModulated, modulationFunction2, 'composition', 'add', 'modulationRegion', modulationRegion);

% oiSequence object for computing a sequence of ois where the oiModulated
% (a grating) is ADDED with the oiBackground using a biphasic modulation function
theOIsequence(3) = oiSequence(oiBackground, oiModulatedGabor, modulationFunction1,  'composition', 'add');

% oiSequence object for computing a sequence of ois where the oiModulated
% (a grating) is BLENDED with the oiBackground using a biphasic modulation function
theOIsequence(4) = oiSequence(oiBackground, oiModulatedGabor, modulationFunction3, 'composition', 'blend');

% Plot the oisequences
colsNum = ceil(theOIsequence(1).length/2);
subplotPosVectors = NicePlot.getSubPlotPosVectors(...
           'rowsNum', 2*numel(theOIsequence), ...
           'colsNum', colsNum+1, ...
           'heightMargin',   0.025, ...
           'widthMargin',    0.02, ...
           'leftMargin',     0.03, ...
           'rightMargin',    0.00, ...
           'bottomMargin',   0.03, ...
           'topMargin',      0.03);
       
hFig = figure(1); clf;
set(hFig, 'Color', [1 1 1], 'Position', [10 10 1350 900]);

for oiIndex = 1:theOIsequence(1).length
    for k = 1:numel(theOIsequence)
        
        row = (k-1)*2+1 + floor((oiIndex-1)/colsNum);
        
        if (oiIndex == 1)
            % Plot the modulation function
            subplot('Position', subplotPosVectors((k-1)*2+1,1).v);
            plot(1:theOIsequence(k).length, theOIsequence(k).modulationFunction, 'rs-', 'LineWidth', 1.5);
            set(gca, 'XLim', [1 theOIsequence(k).length], 'FontSize', 12);
            title(sprintf('composition\n''%s''', theOIsequence(k).composition));
            xlabel('frame index');
            ylabel('modulation');
        end
        
        % Ask theOIsequence to return the oiIndex-th frame
        currentOI = theOIsequence(k).frameAtIndex(oiIndex);
        % And plot it
        col = 2+mod(oiIndex-1, colsNum);
        subplot('Position', subplotPosVectors(row,col).v);
        rgbImage = xyz2srgb(oiGet(currentOI, 'xyz'));
        imagesc(rgbImage, [0 1]);
        title(sprintf('frame %d', oiIndex));
        axis 'image'
        set(gca, 'XTick', [], 'YTick', [],'FontSize', 12);
        drawnow
    end
end
