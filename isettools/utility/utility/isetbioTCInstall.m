function tabComplete = isetbioTCInstall()
% Helps set up auto completion for isetbio objects and functions
%
% Syntax:
%   tabComplete = isetbioTCInstall()
%
% Description:
%	 The settings take effects after reboot the current Matlab session
%
%    The following functions are currently supported:
%       DISPLAY: displayCreate, displayGet, displaySet
%       SCENE:   sceneCreate, sceneSet, sceneGet
%       OI:      oiCreate, oiSet, oiGet
%
% Inputs:
%    None.
%
% Outputs:
%    tabComplete - The XML file containing the created object
%
%  Notes:
%    * [Note: XXX - The basic idea here is to update the TC.xml file.
%      Because we do not have a standard for comments in functions, it's
%      hard to extract potential string options by parsing the function
%      file. The content of auto-completion is hard coded here and if we
%      update those functions some time later, we need to update this file
%      as well.]
%
%  See also:
%    isetbioACUninstall
%

% History:
%    xx/xx/14  HJ   ISETBIO TEAM, 2014
%    12/11/16  dhb  Remove SENSOR options. Think I did this right, but was
%                   not sure how to test.
%    12/14/17  jnm  Formatting

%% Init and backup
pathTC = fullfile(matlabroot, 'toolbox', 'local', 'TC.xml');
backTC = fullfile(matlabroot, 'toolbox', 'local', 'TC_backup.xml');

% back up
if exist(backTC, 'file'), error('backup file already exist'); end
if ~copyfile(pathTC, backTC, 'f'), error('Cannot backup tc file'); end
try fileattrib(pathTC, '+w'); catch, end

tabComplete = xmlread(pathTC);

%% Tab-completion for DISPLAY
%  displayCreate
val  = 'CRT-Dell CRT-HP Dell-Chevron LCD-Apple LCD-Dell LCD-HP OLED-Sony';
node = xmlElementGet(tabComplete, 'displayCreate');
nodeArgValueSet(node, 1, 'VAR', val, tabComplete);

%  displayGet
val  = ['type name gammaTable dacsize nlevels levels ' ...
        'wave nwave spd whiteSPD nprimaries rgb2xyz rgb2lms ' ...
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
       'macbethtungsten macbethEE_IR reflectanceChart ringsRays ' ...
       'harmonic sweepFrequency lineD65 lineEE barEE pointArray ' ...
       'gridLines checkerboard frequencyOrientation slantedEdge ' ...
       'moireOrient zonePlate starPattern letter font whitenoise ' ...
       'linearIntensityRamp uniformEqualEnergy uniformEqualPhoton ' ...
       'uniformBB vernier'];
node = xmlElementGet(tabComplete, 'sceneCreate');
nodeArgValueSet(node, 1, 'VAR', val, tabComplete);

% sceneGet
val = ['name type filename rows cols size height width diagonalsize ' ...
       'heightAndWidth area distance fov hfov vfov aspectratio ' ...
       'magnification depthmap dangular photons knownReflectance ' ...
       'peakRadiance peakRadianceAndWave dataMax dataMin energy ' ...
       'meanEnergySPD meanPhotonsSPD roiPhotonsSPD meanLuminance ' ...
       'luminance xyz lms sampleSize spatialResolution sampleSpacing ' ...
       'distancePerDegree DegreesPerDistance DegreesPerSample ' ...
       'spatialSupport angularResolution frequencyResolution ' ...
       'maxFrequencyResolution frequencySupport fSupportX fSupportY ' ...
       'binWidth wave nwave illuminant illuminantName ' ...
       'illuminantEnergy illuminantPhotons illuminantXYZ ' ...
       'illuminantWave illuminantFormat illuminantComment rgbImage'];
node = xmlElementGet(tabComplete, 'sceneGet');

nodeArgValueSet(node, 1, 'VAR', [], tabComplete);
nodeArgValueSet(node, 2, 'VAR', val, tabComplete);

% sceneSet
val = ['name type distance hfov magnification photons wave depthMap ' ...
       'illuminant peakPhotonRadiance illuminantName illuminantEnergy ' ...
       'illuminantPhotons illuminantComment knownReflectance ' ...
       'luminance meanluminance consistency'];
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
val = ['name type filename consistency rows cols size imageDistance ' ...
       'hfov vfov aspectratio height width diagonal heightAndWidth ' ...
       'area centerPixel photons photonsNoise dataMax dataMin ' ...
       'bitDepth energy energyNoise meanIlluminance illuminance xyz ' ...
       'binWidth wave nwave hSpatialResolution wSpatialResolution ' ...
       'sampleSpacing distancePerSample distancePerDegree ' ...
       'spatialSamplingPositions hAngularResolution ' ...
       'wAngularResolution angularResolution frequencySupport ' ...
       'maxFrequencyResolution fsupportX fsupportY depthMap optics ' ...
       'rgbImage'];
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

%% Write to TC file
xmlwrite(pathTC, tabComplete);

end

%% Get xml node by attribute 'name'
function node = xmlElementGet(obj, name)
% Get xml node by attribute 'name'
%
% Syntax:
%   node = xmlElementGet(obj, name)
%
% Description:
%    Retrieve an xml node by its 'name' attribute
%
% Inputs:
%    obj  - The XML object
%    name - The node's name 
%
% Outputs:
%    node - The xml node

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
function node = nodeArgValueSet(node, index, ctype, val, obj)
% Set argument values for one element/function
%
% Syntax:
%   node = nodeArgValueSet(node, index, ctype, val, obj)
%
% Description:
%    Set the argument values for a single element or function. Please
%    remember that in MATLAB index values start at 1.
%
% Inputs
%    node  - The node to modify
%    index - The child node index
%    ctype - The character type of the value
%    val   - The value to assign to the node
%    obj   - The xml object containing the node
%
% Outputs:
%    node  - The modified xml node

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