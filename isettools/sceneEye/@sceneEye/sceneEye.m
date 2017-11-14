classdef sceneEye < hiddenHandle % TL: What does this hiddenHandle mean? I seem to need it to avoid errors. 
    
    %%sceneEye Create a sceneEye object
    %
    %   sceneEye will contain all the information needed to construct a
    %   new PBRT file which we can then render to get an optical
    %   image. Ideally, sceneEye should be analogous to the "scene"
    %   object in ISET and have similar commands. Unlike the ISET
    %   "scene", however, we follow ISETBIO convention and have it as
    %   a MATLAB class.
    %
    %   This code is starting out in service of the eyeModel always
    %   uses the eye model file. It may generalize later.
    %
    %   Workflow
    %     thisScene = sceneEye('name',xxx,'pbrtFile',xxx);
    %     thisScene.accommodation = double
    %     ...
    %     oi = thisScene.render(varargin);
    %
    % See also:
    %
    % Depends on: pbrt2ISET, ISETBIO
    %
    % TL ISETBIO Team, 2017
    
    % TODO: Is there a way to check inputs? For example, eyePos is not a
    % dependent variable, but is instead read in from the PBRT file.
    % However, say the user wants to change the value in their script so
    % they write:
    %
    % myScene = sceneEye('pbrtFile',xxx);
    % myScene.eyePos = [x y z]; 
    %
    % Is there a way to ensure they put in a 3x1 vector for eyePos, other
    % than just rigourous error checking in the code?
    
    properties (GetAccess=public, SetAccess=public)
        
        %NAME The name of the render
        name;
        
        %RES resolution of render (pixels)
        %   Instead of rows/cols we use a general resolution variable. This
        %   is because the eye model can only take equal rows and columns
        %   and the rendered image is always square.
        resolution;
        
        %FOV Field of view of the render in degrees
        %   This value is calculated from the retina distance and the
        %   retina size. This is only a close approximation since the
        %   retina is very slightly curved.
        fov;
        
        %ACCOMMODATION Diopters of accommodation for 550 nm light
        %   We change the properties of the lens to match the desired
        %   accommodation. For example, if we set this to 5 diopters, 550
        %   nm rays from 0.2 meters will be in focus on the retina.
        accommodation;
        
        %PUPILDIAMETER Diameter of the pupil (mm)
        pupilDiameter;
        
        %RETINARADIUS The curvature of the retina in mm
        %   If one imagines the retina as asection of a sphere, this radius
        %   value determines the distance from the edge of the sphere to
        %   its center.  We will not change this most of the time, but
        %   sometimes it is helpful to make the retina very flat in order
        %   to measure certain properties of the eye.
        retinaRadius;
        
        %retinaDistance Distance between the back lens and the retina
        %   We will not change this most of the time, but sometimes it is
        %   helpful to move the retina back and forth, like a camera
        %   sensor, to see things affects like chromatic aberration.
        retinaDistance;
        
        %NUMRAYS Number of rays to shoot per pixel.
        %   This determines the quality of the render and affects the time
        %   spent rendering. This should be a factor of 2. Low quality is
        %   typically 64 or256 rays, high quality is typically 2048 or 4096
        %   rays.
        numRays;
        
        %NUMBOUNCES Number of bounces before ray terminates
        %   This also determines how accurately light is modeled in the
        %   rendering. The amount needed is highly scene dependent.
        %   Typically set to 1 for simple, diffuse scenes. A high value
        %   would be 4-8 for scenes with lots of reflections, caustics, or
        %   glassy materials.
        numBounces;
        
        %NUMCABANDS Number of wavelength samples to take when modeling CA
        %   We shoot extra rays of different wavelengths in order to model
        %   chromatic aberration through the lens system. This determines
        %   the number of samples we take. For example, if we set this to 4
        %   we shoot rays at...
        numCABands;
        
        %EYEPOS Position of the eye within the scene
        %   [x y z]
        eyePos;
        
        %EYETO Point where the eye is looking at
        %   [x y z], the difference between eyeTo and eyePos is the
        %   direction vector that the optical axis is aligned with.
        eyeTo;
        
        %EYEUP Up vector used when building the LookAt transform
        %   [x y z], this is typically [0 0 1] but it depends on how the
        %   eye is oriented. For example, if this was [0 0 -1] the eye
        %   would be "upside down." Some values are not valid, for example
        %   if the eye is looking down the z-axis (eyePos = [0 0 0], eyeTo
        %   = [0 0 1]) then the up vector cannot be [0 0 1].
        eyeUp;
        
        %DEBUGMODE Toggle debug mode.
        %   For debug mode we switch to a perspective camera with the same
        %   FOV as the eye. This can be potentially faster and easier to
        %   render than going through the eye.
        debugMode;
    end
    
    properties (Dependent)
        
        %WIDTH Width of imaged retina (mm)
        %   Depends on fov, retinaDistance, and rows/cols
        width;
        
        %HEIGHT Height of imaged retina (mm)
        %   Depends on fov, retinaDistance, and rows/cols
        height;
        
        %SAMPLESIZE Samples spacing, e.g. width/xRes and height/yRes.
        %   We assume square samples. This is not always accurate at large
        %   fov's.
        sampleSize;
        
    end
    
    properties(GetAccess=public, SetAccess=private)
        
        % Actually...should we include this? This will change depending on
        % the accommodation...
        %LENSFILE Path to the .dat file that describes the lens system
        %   This file includes descriptions of the thickness, curvature,
        %   and diameter of the various components in the eye.
        lensFile;
        
        % PBRTFILE Path to the original .pbrt file this scene is based on
        %   Depends on the pbrt file used to create the scene. Should not
        %   be changed.
        pbrtFile;
        
        %WORKINGDIR Directory to store temp files needed for rendering
        %   We make a copy of the scene into the working directory, and
        %   then output new PBRT files into this directory. We also save
        %   the raw rendered data (xxx.dat) in this folder.
        workingDir;
        
    end
    
    properties(Access=private)
        
        % Pretty much everything else we read in the PBRT file that we
        % don't want them to see but still need in order to write out the
        % new PBRT file. This includes things like the WorldBegin/WorldEnd
        % block, the PixelFilter, maybe the SurfaceIntegrator?
        
        %RECIPE Structure that holds all instructions needed to
        %   render the PBRT file. 
        recipe;
        
    end
    
    properties (Constant)
        %WAVE
        wave = []; % TODO Look it up and fill it in
        
        
    end
    
    methods
        % Constructor
        function obj = sceneEye(pbrtFile,varargin)
            % Initialize the sceneEye class
            %
            % Reads a PBRT file and fills out the information needed
            % for the sceneEye object.  That object will be rendered
            % using the PBRT methods (docker image).
            
            p = inputParser;
            p.KeepUnmatched = true;
            p.addRequired('pbrtFile',@ischar); % Either a pbrt file or just a scene name
            p.addParameter('name','scene-001',@ischar);
            p.addParameter('workingDirectory','',@ischar);
             
            % An optional parameter used by scenes that consist of only a
            % planar surface (e.g. slanted bar). We will move the plane to
            % the given distance in mm.
            p.addParameter('planeDistance',1000,@isnumeric);
            
            p.parse(pbrtFile,varargin{:});
            
            % Read in PBRT file
            [~,name,ext] = fileparts(pbrtFile);
            
            if(isempty(ext))
                % The user has given us a scene name and not a full pbrt
                % file. Let's find the right file.
                switch name
                    case('numbersAtDepth')
                        scenePath = fullfile(isetbioDataPath,'pbrtscenes', ...
                            'NumbersAtDepth','numbersAtDepth.pbrt');
                    case('slantedBar')
                        scenePath = fullfile(isetbioDataPath,'pbrtscenes', ...
                            'SlantedBar','slantedBar.pbrt');
                    case('chessSet')
                        scenePath = fullfile(isetbioDataPath,'pbrtscenes', ...
                            'ChessSet','chessSet.pbrt');
                    otherwise
                        error('Did not recognize scene type.');
                end
            else
                scenePath = pbrtFile;
            end
            
            % Setup working folder
            if(isempty(p.Results.workingDirectory))
                [path,~,~] = fileparts(scenePath); % Determine scene folder name from scene path
                [~,sceneFolder] = fileparts(path);
                obj.workingDir = fullfile(isetbioRootPath,'local',sceneFolder);
            else
                obj.workingDir = p.Results.workingDirectory;
            end
            
            obj.pbrtFile = createWorkingFolder(scenePath,'workingDir',obj.workingDir);
            
            % Parse PBRT file
            recipe = piRead(obj.pbrtFile);
            recipe.outputFile = obj.pbrtFile;
            
            % Move planar scenes (e.g. slantedBar) to the desired distance
            if(strcmp(name,'slantedBar'))
                recipe = piMoveObject(recipe,'1_WhiteCube','Translate',[0 p.Results.planeDistance 0]);
                recipe = piMoveObject(recipe,'2_BlackCube','Translate',[0 p.Results.planeDistance 0]);
            end
                                    
            % Note: What happens if the recipe doesn't include any of
            % the following, or any of the subfields we call?
            
            % Check to make sure this PBRT file has a realistic eye.
            if(~strcmp(recipe.camera.subtype,'realisticEye'))
                error('This PBRT file does not include a "realistic eye" camera class.')
                % TODO: Overwrite the camera class with a realistic eye if
                % this happens.
            end
            
            % Set properties
            obj.name = p.Results.name;
            obj.resolution = recipe.film.xresolution.value;
            obj.retinaDistance = recipe.camera.retinaDistance.value;
            obj.pupilDiameter = recipe.camera.pupilDiameter.value;
             
            obj.retinaDistance = recipe.camera.retinaDistance.value;
            obj.retinaRadius = recipe.camera.retinaRadius.value;
            
            retinaSemiDiam = recipe.camera.retinaSemiDiam.value;
            obj.fov = 2*atand(retinaSemiDiam/obj.retinaDistance);
            
            % There's no variable for accommodation but we can infer it
            % from the name of the lens. We assume the naming conventions
            % is "%s_%f.dat" This is not foolproof, so maybe we can think
            % of a more robust way to do this in the future?
            obj.lensFile = recipe.camera.specfile.value;
            % Use regular expressions to find any floats within the string
            value = regexp(obj.lensFile, '(\d+,)*\d+(\.\d*)?', 'match');
            obj.accommodation = str2double(value{1});
            
            obj.numRays = recipe.sampler.pixelsamples.value;
            
            % These two are often empty, so let's do checks here. However,
            % I should find a more permanant solution to cases like these.
            % (See note above).
            % Maybe in piGetRenderRecipe we should put in the default
            % values if any of these rendering options are missing
            % (e.g. if Renderer is missing, put in Renderer 'sampler'.)
            if(isfield(recipe.integrator,'maxdepth'))
                obj.numBounces = recipe.integrator.maxdepth.value;
            else
                obj.numBounces = 1;
            end
            if(isfield(recipe.renderer,'nWaveBands'))
                obj.numCABands = recipe.renderer.nWaveBands.value;
            else
                obj.numCABands = 0;
            end
            
            % What happens if there is no look at, only a Transform?
            % TODO:
            % 1. Read Look At
            % 2. Transform it to find new Look At
            % 3. Populate eyePos/eyeTo/eyeUp using Look At
            % else
            % 1. Directly populate eyePos/eyeTo/eyeUp using Look At
            if(~isempty(recipe.lookAt))
                obj.eyePos = recipe.lookAt.from;
                obj.eyeTo = recipe.lookAt.to;
                obj.eyeUp = recipe.lookAt.up;
            end
            
            obj.recipe = recipe;
            obj.debugMode = false;
            
        end
        
        %% Get methods for dependent variables
        
        function val = get.width(obj)
            % Rendered image is alway square.
            val = 2*tand(obj.fov/2)*obj.retinaDistance;
        end
        
        function val = get.height(obj)
            % Rendered image is alway square.
            val = 2*tand(obj.fov/2)*obj.retinaDistance;
        end
        
        function val = get.sampleSize(obj)
            val = (2*tand(obj.fov/2)*obj.retinaDistance)/obj.resolution;
        end
        
        %% Set methods for dependent variables
        
        function set.width(obj, val)
            obj.fov = 2*atand((val/2)/obj.retinaDistance);
        end
        
        function set.height(obj, val)
            obj.fov = 2*atand((val/2)/obj.retinaDistance);
        end
        
        function set.sampleSize(obj, val)
            obj.fov = 2*atand((val*obj.resolution/2)/obj.retinaDistance);
        end
        
        
    end
    
    methods (Access=public)
        
        [oi, terminalOutput, outputFile] = render(obj,varargin);
        
        
    end
        
   
    
end

