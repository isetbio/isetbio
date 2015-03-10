function tabComplete = isetbioTCInstall()
%% function tabComplete = isetbioTCInstall
%    This function helps set up auto completion for isetbio objects and
%    functions
%
%    The settings take effects after reboot the current Matlab session
%
%  Functions currently supported:
%    DISPLAY: displayCreate, displayGet, displaySet
%    SCENE:   sceneCreate, sceneSet, sceneGet
%    OI:      oiCreate, oiSet, oiGet
%    SENSOR:  sensorCreate, sensorGet, sensorSet
%  
%  Might add for future:
%    PIXEL:   pixelSet, pixelGet
%    MACULAR: macularGet, macularSet
%    LENS:    lensGet, lensSet
%    EYEMOVE: emGet, emSet
%  
%  Notes:
%    The basic idea here is to update the TC.xml file. Because we do not
%    have a standard for comments in functions, it's hard to extract
%    potential string options by parsing the function file. The content of
%    auto-completion is hard coded here and if we update those functions
%    some time later, we need to update this file as well.
%
%  See also:
%    isetbioACUninstall
%
%  (HJ) ISETBIO TEAM, 2014

%% Init and backup
pathTC = fullfile(matlabroot, 'toolbox', 'local', 'TC.xml');
backTC = fullfile(matlabroot, 'toolbox', 'local', 'TC_backup.xml');

% back up
if exist(backTC, 'file'), error('backup file already exist'); end
if ~copyfile(pathTC, backTC, 'f'), error('Cannot backup tc file'); end
try fileattrib(pathTC,'+w'); catch, end

tabComplete = xmlread(pathTC);

%% Tab-completion for DISPLAY
%  displayCreate
val  = 'CRT-Dell CRT-HP Dell-Chevron LCD-Apple LCD-Dell LCD-HP OLED-Sony';
node = xmlElementGet(tabComplete, 'displayCreate');
nodeArgValueSet(node, 1, 'VAR', val, tabComplete);

%  displayGet
val  = ['type name gammaTable dacsize nlevels levels ' ...
        'wave nwave spd whiteSPD nprimaries rgb2xyz rgb2lms' ...
        'whiteXYZ whitexy whiteLMS primariesXYZ priamriesxy ' ...
        'dpi ppi metersPerDot dotsPerMeter dotsPerDeg viewingDistance'];
node = xmlElementGet(tabComplete, 'displayGet');

nodeArgValueSet(node, 1, 'VAR', [], tabComplete);
nodeArgValueSet(node, 2, 'VAR', val, tabComplete);

% displaySet
val = 'psf gTable wave spd dpi ppi psf viewingDistance comment';
node = xmlElementGet(tabComplete, 'displaySet');

nodeArgValueSet(node, 1, 'VAR', [], tabComplete);
nodeArgValueSet(node, 2, 'VAR', val, tabComplete);
nodeArgValueSet(node, 3, 'VAR', [], tabComplete);

%% Tab-Completion for SCENE
%  sceneCreate
val = ['default macbethd65 macbethd50 macbethillc macbethfluorescent ' ...
       'macbethtungsten macbethEE_IR reflectanceChart ringsRays '...
       'harmonic sweepFrequency lineD65 lineEE barEE pointArray '...
       'gridLines checkerboard frequencyOrientation slantedEdge ' ...
       'moireOrient zonePlate starPattern letter font whitenoise ' ...
       'linearIntensityRamp uniformEqualEnergy uniformEqualPhoton ' ...
       'uniformBB vernier'];
node = xmlElementGet(tabComplete, 'sceneCreate');
nodeArgValueSet(node, 1, 'VAR', val, tabComplete);

% sceneGet
val = ['name type filename rows cols size height width diagonalsize '...
       'heightAndWidth area distance fov hfov vfov aspectratio '...
       'magnification depthmap dangular photons knownReflectance '...
       'peakRadiance peakRadianceAndWave dataMax dataMin energy '...
       'meanEnergySPD meanPhotonsSPD roiPhotonsSPD meanLuminance '...
       'luminance xyz lms sampleSize spatialResolution sampleSpacing '...
       'distancePerDegree DegreesPerDistance DegreesPerSample '...
       'spatialSupport angularResolution frequencyResolution '...
       'maxFrequencyResolution frequencySupport fSupportX fSupportY '...
       'binWidth wave nwave illuminant illuminantName illuminantEnergy '...
       'illuminantPhotons illuminantXYZ illuminantWave illuminantFormat'...
       'illuminantComment rgbImage'];
node = xmlElementGet(tabComplete, 'sceneGet');

nodeArgValueSet(node, 1, 'VAR', [], tabComplete);
nodeArgValueSet(node, 2, 'VAR', val, tabComplete);

% sceneSet
val = ['name type distance hfov magnification photons wave depthMap '...
       'illuminant peakPhotonRadiance illuminantName illuminantEnergy '...
       'illuminantPhotons illuminantComment knownReflectance luminance'...
       'meanluminance consistency'];
node = xmlElementGet(tabComplete, 'sceneSet');

nodeArgValueSet(node, 1, 'VAR', [], tabComplete);
nodeArgValueSet(node, 2, 'VAR', val, tabComplete);
nodeArgValueSet(node, 3, 'VAR', [], tabComplete);

%% Tab-Completion for Optical Image
%  oiCreate
val = 'default uniformd65 uniformee human wvfHuman';
node = xmlElementGet(tabComplete, 'oiCreate');

nodeArgValueSet(node, 1, 'VAR', val, tabComplete);

% oiGet
val = ['name type filename consistency rows cols size imageDistance '...
       'hfov vfov aspectratio height width diagonal heightAndWidth '...
       'area centerPixel photons photonsNoise dataMax dataMin '...
       'bitDepth energy energyNoise meanIlluminance illuminance xyz '...
       'binWidth wave nwave hSpatialResolution wSpatialResolution '...
       'sampleSpacing distancePerSample distancePerDegree '...
       'spatialSamplingPositions hAngularResolution wAngularResolution '...
       'angularResolution frequencySupport maxFrequencyResolution '...
       'fsupportX fsupportY depthMap optics rgbImage'];
node = xmlElementGet(tabComplete, 'oiGet');

nodeArgValueSet(node, 1, 'VAR', [], tabComplete);
nodeArgValueSet(node, 2, 'VAR', val, tabComplete);

% oiSet
val = ['name distance hfov magnification photons wave optics ' ...
       'consistency gamma'];
node = xmlElementGet(tabComplete, 'oiSet');

nodeArgValueSet(node, 1, 'VAR', [], tabComplete);
nodeArgValueSet(node, 2, 'VAR', val, tabComplete);
nodeArgValueSet(node, 3, 'VAR', [], tabComplete);

%% Tab-Completion for SENSOR
%  sensorCreate
val = ['bayer-grbg bayer-rggb bayer-bggr bayer-gbrg bayer-ycmy ' ...
       'bayer-cyym ideal monochrome grbc interleaved fourColor custom' ...
       'human'];
node = xmlElementGet(tabComplete, 'sensorCreate');

nodeArgValueSet(node, 1, 'VAR', val, tabComplete);

% sensorGet
val = ['name type row col size height width dimension spatialSupport '...
       'wSpatialSupport hSpatialSupport hfov vfov hDegPerPixel '...
       'vDegPerPixel hDegPerDistance vDegPerDistance fov microLens '...
       'chiefRayAngle chiefRayAngleDegrees sensorEtendue volts '...
       'digitalValues electrons photons roiLocs roiRect roiVolts '...
       'roiElectrons roiVoltsMean roiElectronsMean hLineVolts '...
       'vLineVolts hLineElectrons vLineElectrons responseRatio roi '...
       'analogGain analogOffset sensorDynamicRange quantization wave '...
       'binWidth nWave filterTransmissivities infraredFilter cfaName '...
       'filterNames nfilters filterColorLetters filterPlotColors '...
       'filterColorLettersCell spectralQE pattern dsnuSigma prnuSigma '...
       'fpnParameters dsnuImage prnuImage columnfpn columnDsnu '...
       'columnPRNU colOffsetFpnVector colGainFpnVector noiseFlag '...
       'reuseNoise noiseSeed pixel autoExposure exposureTime '...
       'uniqueExpTimes exposurePlane cds nSamplesPerPixel consistency '...
       'human humanLens humanLensTransmittance humanLensAbsorption '...
       'humanMacular humanMacularDensities humanMacularTransimittance '...
       'humanMacularAbsorption humanOcularTransmittance coneType '...
       'humanConeDensity xy adaptationGain adaptationOffset '...
       'addaptedVolts timeInterval emType emTremor emDrift totalTime '...
       'emMicrosaccade rgb'];
node = xmlElementGet(tabComplete, 'sensorGet');

nodeArgValueSet(node, 1, 'VAR', [], tabComplete);
nodeArgValueSet(node, 2, 'VAR', val, tabComplete);

% sensorSet
val = ['name rows cols size fov filterSpectra filterNames wave pattern '...
       'infaredFilter patternAndSize volts digitalValues analogGain '...
       'analogOffset roi cds quantizationMethod dsnuImage prnuImage '...
       'dsnuLevel prnuLevel columnFpnParameters colGainFpnVector '...
       'colOffsetFpnVector noiseFlag reuseNoise noiseSeed exposureTime '...
       'exposurePlane autoExposure pixel pixelVignetting microLens '...
       'sensorEtendue microLensOffset sensorComputeMethod ngridSamples '...
       'consistency mccRectHandles mccCornerPoints sensorPositions '...
       'humanLens humanCone humanConeType humanConeDensities totalTime '...
       'adaptationGain adaptationOffset timeInterval emType emTremor '...
       'emDrift emMicrosaccade'];
node = xmlElementGet(tabComplete, 'sensorSet');

nodeArgValueSet(node, 1, 'VAR', [], tabComplete);
nodeArgValueSet(node, 2, 'VAR', val, tabComplete);
nodeArgValueSet(node, 3, 'VAR', [], tabComplete);

%% Write to TC file
xmlwrite(pathTC, tabComplete);

end

%% Get xml node by attribute 'name'
function node = xmlElementGet(obj, name)
    allItems = obj.getElementsByTagName('binding');
    for ii = 0 : allItems.getLength - 1
        if strcmp(allItems.item(ii).getAttribute('name'), name)
            warning('settings for %s already exists', name);
            node = allItems.item(ii);
            return;
        end
    end
    node = obj.createElement('binding');
    node.setAttribute('name', name);
    obj.getLastChild.appendChild(node);
end

%% Set argument values for one element / function
%  index is starting from 1
function node = nodeArgValueSet(node, index, ctype, val, obj)
    arg = node.getChildNodes.item(index - 1); % get argument
    if isempty(arg)
        arg = obj.createElement('arg');
        node.appendChild(arg);
    end
    arg.setAttribute('argn', num2str(index));
    arg.setAttribute('ctype', ctype);
    if ~isempty(val)
        arg.setAttribute('value', val);
    end
end