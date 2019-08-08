function [recipe, sceneUnits, workingDir, origPath] = ...
    loadPbrtScene(pbrtFile, varargin)
% Set up PBRT scene given a name or file location. See description below.
%
% Syntax:
%   [recipe, sceneUnits, workingDir, origPath] = ...
%       loadPbrtScene(pbrt, [varargin])
%
% Description:
%    Setup a PBRT scene given it's name or file location. Primarily
%    includes the following steps:
%      1. Check if we're given a pbrt file or a scene name.
%         a) If a pbrt file, read it, and return a recipe
%         b) If a scene name, download it from the RDT, read it, and return
%            a recipe.
%      2. Set up a working folder derived from the scene name. Copy all
%         necessary files over to the newly created working directory.
%      3. Apply any adjustable parameters given by the user to the recipe,
%         e.g. moving a planar target a certain distance away.
%
%    The user can call sceneEye with the name of a scene to automatically
%    load it. To remove bloat from teh actual sceneEye class, we will do
%    that parsing/selection in this function instead.
%
% Inputs:
%    sceneName   - String. Either a scene name like "slantedBar" or an
%                  actual pbrt filepath like ("xxx.pbrt")
%    p           - Object. An inputParser from sceneEye needs to be passed
%                  in so we can find certain parameters (e.g.
%                  planeDistance) when setting up the scene.
%
% Outputs:
%    recipe      - Object. A recipe object of selected scene after
%                  adjustment parameters are applied.
%    sceneUnits  - String. some scenes are in meters and some are in
%                  millimeters. There is a flag in the sceneEye class to
%                  specify this and render appropriately.
%    workingDir  - String. A created working directory derived from the
%                  scene name.
%    origPath    - String. The original path to the pbrt file.
%
% Optional key/value pairs:
%    name        - String. A string containing the desired scene name.
%                  Default 'scene-001'.
%    workingDirectory
%                - String. A string containing the original working
%                  directory. Default ''.
%    planeDistance
%                - Numeric. The distance to the plane. Default 1.
%    planeTexture
%                - String. A string filepath to a texture file for the
%                  plane. Default is squareResolutionChart.exr inside the
%                  data & imageTextures directory off the root path.
%    planeSize   - Matrix. A 1x2 matrix of the height and width of the
%                  plane. Default is [1 1].
%    pointDiameter
%                - Numeric. The diameter of a point. Default is 0.001.
%    pointDistance
%                - Numeric. The distance to the point. Default 1.
%    gamma       - Boolean. Whether to use sRGB or 1RGB. Default true(sRGB)
%    useDisplaySPD
%                - Numeric. A numeric boolean indicating whether or not to
%                  use the display SPD. Default 0 (no).
%    topDepth    - Numeric. A slanted bar texture parameter of the top
%                  depth. Default 1.
%    bottomDepth - Numeric. A slanted bar texture parameter of the bottom
%                  depth. Default 1.
%    eccentricity
%                - Numeric. A slanted bar parameter representing the
%                  eccentricity. Default 0.
%    objectDistance
%                - Numeric. A snellen single parameter for the object
%                  distance. Default 1.
%    objectSize  - Numeric. A snellen single parameter for the object size.
%                  Default 0.3.
%    Adist       - Numeric. A letters at depth parameter, representing the
%                  distance to the 'A'. Default 0.1.
%    Bdist       - Numeric. A letters at depth parameter, representing the
%                  distance to the 'B'. Default 0.2.
%    Cdist       - Numeric. A letters at depth parameter, representing the
%                  distance to the 'C'. Default 0.3.
%    Adeg        - Numeric. A letters at depth parameter, representing the
%                  degrees of distance from the middle for the letter A.
%                  Default 6.
%    Cdeg        - Numeric. A letters at depth parameter, representing the
%                  degrees of distance from the middle for the letter C.
%                  Default 4.
%    nChecks     - Matrix. A 1x2 matrix. This is a letters at depth
%                  parameter, representing the number of checks along the
%                  wall and ground in the background. Default [64 64].
%
% Notes:
%    * TODO: I'd like to keep splitting up the above steps into more
%      functions to neaten things up. I think we might even be able to
%      combine step 1 and step 3.
%

% History:
%    05/25/18  TL   Created
%    05/13/19  JNM  Documentation pass
%    08/08/19  JNM  Merge master in

%% Parse inputs
% Because the inputs have been passed in first through sceneEye, and then
% through loadPbrtScene, they seem to be in a nested cell. 
varargin = varargin{:};

p = inputParser;
p.addRequired('pbrtFile', @ischar);

p.addParameter('name', 'scene-001', @ischar);
p.addParameter('workingDirectory', '', @ischar);

% The following are optional parameters used by unique scenes (e.g.
% slantedBar, texturedPlane, pointSource). We can use these parameters to
% move the plane/point to the given distance (in mm) and, if applicable,
% attach the provided texture.
p.addParameter('planeDistance', 1, @isnumeric);
p.addParameter('planeTexture', fullfile(piRootPath, 'data', ...
    'imageTextures', 'squareResolutionChart.exr'), @ischar);
p.addParameter('planeSize', [1 1], @isnumeric);
p.addParameter('pointDiameter', 0.001, @isnumeric);
p.addParameter('pointDistance', 1, @isnumeric);

% texturedPlane
p.addParameter('gamma', 'true', @ischar);
p.addParameter('useDisplaySPD', 0);

% slantedBarTexture
p.addParameter('topDepth', 1, @isnumeric); 
p.addParameter('bottomDepth', 1, @isnumeric); 

% slantedBar
p.addParameter('eccentricity', 0, @isnumeric); 

% snellenSingle
p.addParameter('objectDistance', 1, @isnumeric); 
p.addParameter('objectSize', 0.3, @isnumeric); 

%lettersAtDepth and lettersAtDepthPlus
p.addParameter('Adist', 0.1, @isnumeric); 
p.addParameter('Bdist', 0.2, @isnumeric); 
p.addParameter('Cdist', 0.3, @isnumeric); 
p.addParameter('Ddist', [], @isnumeric); 
p.addParameter('Edist', [], @isnumeric); 
p.addParameter('Adeg', 6, @isnumeric); 
p.addParameter('Bdeg', [], @isnumeric);
p.addParameter('Cdeg', 4, @isnumeric); 
p.addParameter('Ddeg', [], @isnumeric);
p.addParameter('Edeg', [], @isnumeric);
p.addParameter('nchecks', [64 64], @isnumeric); % [wall, ground]

p.parse(pbrtFile, varargin{:});

% Default
sceneUnits = 'm';

%% Check if we've been given a sceneName or a pbrt file.
[~, sceneName, ext] = fileparts(pbrtFile);
if(isempty(ext))
    % scene name
    sceneNameFlag = true;
else
    % pbrt file
    sceneNameFlag = false;
    scenePath = pbrtFile;
end

%% Load the scene
if(sceneNameFlag)
    % The user has given us a scene name and not a full pbrt file. Let's
    % find or download the right file.
    switch sceneName
        case('colorfulScene')
            scenePath = fullfile(piRootPath, 'local', 'scenes', ...
                'colorfulScene', 'ColorfulScene.pbrt');
            sceneUnits = 'm';
            pullSceneFromRDT('ColorfulScene', scenePath);

        case('figure')
            scenePath = fullfile(piRootPath, 'data', ...
                'V3', 'figure', 'figure.pbrt');
            sceneUnits = 'm';

        case('coloredCube')
            scenePath = fullfile(piRootPath, 'data', ...
                'V3', 'coloredCube', 'coloredCube.pbrt');
            sceneUnits = 'm';

        case('snellenSingle')
            scenePath = fullfile(piRootPath, 'local', 'scenes', ...
                'snellenSingle', 'snellen_single.pbrt');
            sceneUnits = 'm';
            pullSceneFromRDT('snellenSingle', scenePath);

        case ('snellenAtDepth')
            scenePath = fullfile(piRootPath, 'local', 'scenes', ...
                'snellenAtDepth', 'snellen.pbrt');
            sceneUnits = 'm';
            pullSceneFromRDT('snellenAtDepth', scenePath);

        case ('blackBackdrop')
            scenePath = fullfile(piRootPath, 'data', 'V3', ...
                'blackBackdrop', 'blackBackdrop.pbrt');
            sceneUnits = 'm';

        case ('blankScene')
            scenePath = fullfile(piRootPath, 'data', 'V3', ...
                'blankScene', 'blankScene.pbrt');
            sceneUnits = 'm';

        case('numbersAtDepth')
            scenePath = fullfile(piRootPath, 'local', ...
                'scenes', 'NumbersAtDepth', 'numbersAtDepth.pbrt');
            sceneUnits = 'm';
            pullSceneFromRDT('NumbersAtDepth', scenePath);

        case('slantedBar')
            scenePath = fullfile(piRootPath, 'data', ...
                'V3', 'slantedBar', 'slantedBar.pbrt');
            sceneUnits = 'm';

        case('chessSet')
            scenePath = fullfile(piRootPath, 'local', 'scenes', ...
                'ChessSet', 'chessSet.pbrt');
            sceneUnits = 'm';
            pullSceneFromRDT('ChessSet', scenePath);

        case('chessSet-2')
            % Chess set but with a couple of different materials. This
            % scene can probably go away in the future when it becomes
            % easier to adjust the materials directly within the script.
            scenePath = fullfile(piRootPath, 'local', 'scenes', ...
                'ChessSet_2', 'chessSet2.pbrt');
            sceneUnits = 'm';
            pullSceneFromRDT('ChessSet_2', scenePath);

        case('chessSetScaled')
            scenePath = fullfile(piRootPath, 'local', 'scenes', ...
                'ChessSetScaled', 'chessSetScaled.pbrt');
            sceneUnits = 'm';
            pullSceneFromRDT('ChessSetScaled', scenePath);

        case('texturedPlane')
            scenePath = fullfile(piRootPath, 'data', ...
                'V3', 'texturedPlane', 'texturedPlane.pbrt');
            sceneUnits = 'm';

        case('pointSource')
            scenePath = fullfile(piRootPath, 'data', 'V3', ...
                'SimplePoint', 'simplePointV3.pbrt');
            sceneUnits = 'm';

        case('slantedBarTexture')
            % A variation of slantedBar where planes have texture.
            scenePath = fullfile(piRootPath, 'data', ...
                'V3', 'slantedBarTexture', 'slantedBarTexture.pbrt');
            sceneUnits = 'm';

        case('lettersAtDepth')
            % A, B, C placed at different depths. The depths will be
            % adjustable.
            scenePath = fullfile(piRootPath, 'local', 'scenes', ...
                'lettersAtDepth', 'lettersAtDepth.pbrt');
            sceneUnits = 'm';
            pullSceneFromRDT('lettersAtDepth',scenePath);

        case('lettersAtDepthPlus')
            % A, B, C, D and E at different depths.
            
            scenePath = fullfile(piRootPath, 'local', 'scenes',...
                'lettersAtDepthPlus', 'lettersAtDepthPlus.pbrt');
            sceneUnits = 'm';
            pullSceneFromRDT('lettersAtDepthPlus', scenePath);

        otherwise
            error('Did not recognize scene type.');
    end
end

%% Read the filename and get a recipe
recipe = piRead(scenePath, 'version', 3);
recipe.inputFile = scenePath;

%% Setup the working folder
if(isempty(p.Results.workingDirectory))
    % Determine scene folder name from scene path
    [path, ~, ~] = fileparts(scenePath);
    [~, sceneFolder] = fileparts(path);
    workingDir = fullfile(isetbioRootPath, 'local', sceneFolder);
else
    workingDir = p.Results.workingDirectory;
end

% Copy contents of the working folder over to the local folder.
% We may no longer need to do this, I think it is taken care of in
% recipe.set. Right now it might just be redundant. 
origPath = createWorkingFolder(scenePath, 'workingDir', workingDir);

%% Make adjustments to the recipe
% E.g. move the plane to a certain distance
if(sceneNameFlag)
    switch sceneName
        case('lettersAtDepth')
            % Move the letters in the scene. To do this, we're actually
            % going to remake the scene.
            recipe = piCreateLettersAtDepth('Adist', p.Results.Adist, ...
                'Bdist', p.Results.Bdist, ...
                'Cdist', p.Results.Cdist, ...
                'Adeg', p.Results.Adeg, ...
                'Cdeg', p.Results.Cdeg, ...
                'nchecks', p.Results.nchecks);

        case('lettersAtDepthPlus')
            % Move the letters in the scene. To do this, we're actually
            % going to remake the scene.
            recipe = piCreateLettersAtDepthPlus(...
                'Adist', p.Results.Adist, ...
                'Bdist', p.Results.Bdist, ...
                'Cdist', p.Results.Cdist, ...
                'Ddist', p.Results.Ddist, ...
                'Edist', p.Results.Edist, ...
                'Adeg', p.Results.Adeg, ...
                'Bdeg', p.Results.Bdeg, ...
                'Ddeg', p.Results.Ddeg, ...
                'Edeg', p.Results.Edeg, ...
                'nchecks', p.Results.nchecks);

        case('slantedBar')
            recipe = piCreateSlantedBarScene(...
                'planeDepth', p.Results.planeDistance, ...
                'eccentricity', p.Results.eccentricity);

        case('slantedBarTexture')
            recipe = piCreateSlantedBarTextureScene(...
                'topDepth', p.Results.topDepth, ...
                'bottomDepth', p.Results.bottomDepth);

        case('pointSource')
            % Clear previous transforms
            piClearObjectTransforms(recipe, 'Point');
            piClearObjectTransforms(recipe, 'Plane');
            % Add given transforms
            recipe = piObjectTransform(recipe, 'Point', 'Scale', ...
                [p.Results.pointDiameter p.Results.pointDiameter 1]);
            recipe = piObjectTransform(recipe, 'Point', 'Translate', ...
                [0 0 p.Results.pointDistance]);
            % Make it large!
            recipe = piObjectTransform(recipe, 'Plane', 'Scale', ...
                [p.Results.pointDistance * 10, ...
                p.Results.pointDistance * 10 1]);
            % Move it slightly beyond the point
            recipe = piObjectTransform(recipe, 'Plane', 'Translate', ...
                [0 0 p.Results.pointDistance + 0.5]);

        case('snellenSingle')
            scaling = [p.Results.objectSize p.Results.objectSize 1] ./ ...
                [1 1 1];
            recipe = piObjectTransform(recipe, 'Snellen', ...
                'Scale', scaling);
            recipe = piObjectTransform(recipe, 'Snellen', ...
                'Translate', [0 0 p.Results.objectDistance]);

        case('texturedPlane')
            % Scale and translate
            planeSize = p.Results.planeSize;
            scaling = [planeSize(1) planeSize(2) 1] ./ [1 1 1];
            recipe = piObjectTransform(recipe, 'Plane', 'Scale', scaling);
            recipe = piObjectTransform(recipe, 'Plane', ...
                'Translate', [0 0 p.Results.planeDistance]);
            % Texture
            [pathTex, nameTex, extTex] = fileparts(p.Results.planeTexture);
            copyfile(p.Results.planeTexture, workingDir);
            if(isempty(pathTex))
                error('Image texture must be an absolute path.');
            end
            recipe = piWorldFindAndReplace(recipe, 'dummyTexture.exr', ...
                strcat(nameTex, extTex));

            % If true, use the lcd-apple display primaries to convert to
            % RGB texture values to spectra.
            if(p.Results.useDisplaySPD)
                recipe = piWorldFindAndReplace(recipe, ...
                    '"bool useSPD" "false"', '"bool useSPD" "true"');
            end

            % If true, we convert from sRGB to lRGB in PBRT.
            if(strcmp(p.Results.gamma, 'false'))
                recipe = piWorldFindAndReplace(recipe, ...
                    '"bool gamma" "true"', '"bool gamma" "false"');
            end
    end
end

