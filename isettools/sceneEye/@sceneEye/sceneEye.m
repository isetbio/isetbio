classdef sceneEye < hiddenHandle
% Create a sceneEye object
%
% Syntax:
%   myScene = sceneEye();
%
% Description:
%    sceneEye contains the information needed to construct a new PBRT
%    file that we can then render to get a retinal image.
%
%    sceneEye is analogous to the "scene" structure in ISETBIO (and
%    ISET), and it will support similar commands. Unlike the
%    ISET/ISETBIO "scene", as for new entities we created it as a
%    MATLAB class.
%
%    This code is starting out in service of the eyeModel. We may
%    extend to replace the scene structure (some day).
%
% Notes:
%    * TODO - Implement the check that BW describes below. Is there a way
%      to check inputs? For example, eyePos is not a dependent variable, 
%      but is instead read in from the PBRT file. However, say the user
%      wants to change the value in their script so they write:
%           myScene = sceneEye('pbrtFile', xxx);
%           myScene.eyePos = [x y z];
%      Is there a way to ensure they put in a 3x1 vector for eyePos, other
%      than just rigourous error checking in the code?
%	 * [Note: BW - (Reference first TODO) Yes, I think so, using the
%	   myScene.set('eyePos', val) approach, or perhaps myScene.set.eyePos =
%	   val. In these cases the set operation can pass through an input
%	   parser that validates the input value (I think).]
%    * [Note: BW - maxDepth & nWaveBands are often empty, so let's perform
%      the checks below. However, I should find a more permanant solution
%      to cases like these. (See the note above). Maybe in
%      piGetRenderRecipe we should put in the default values if any of
%      these rendering options are missing (e.g. if Renderer is missing, 
%      put in Renderer 'sampler'.)]
%    * [Note: XXX - (from constructor) What happens if the recipe doesn't
%      include any of the following, or any of the subfields we call?]
%    * TODO - Determine a better way to infer the accommodation. Currently
%      we assume the naming conventions is "%s_%f.dat" This is not
%      foolproof, so maybe we can think of a more robust way to do this in
%      the future?
%    * [Note: TL - What does this hiddenHandle mean? I seem to need it to
%      avoid errors.]
%    * TODO - Fix example!
%
% See Also:
%    Dependencies: pbrt2ISET, ISETBIO
%

% History:
%    xx/xx/17  TL   ISETBIO Team, 2017
%    12/19/17  jnm  Formatting
%    08/08/19  JNM  Merge master in

% Examples:
%{
    % ETTBSkip.  Skip this example in ETTB, since it is known not to work.
    % When the example gets fixed, remove this line and the one above.

    scene3d = sceneEye('chessSet');
               
    scene3d.fov = 30; 
    scene3d.resolution = 128;
    scene3d.numRays = 128;
    scene3d.numCABands = 0;
    scene3d.accommodation = 1; 

    oi = scene3d.render();
    ieAddObject(oi);
    oiWindow;

%}

properties (GetAccess=public, SetAccess=public)
    % name - The name of the render
    name;
    
    % modelname - The name of the schematic eye used to render the scene.
    %   Depending on the model chosen, some other options may not be
    %   applicable. Currently possible models include: navarro (default)
    %   le grand, and arizona eye model.
    modelName;
    
    % resolution - resolution of render (pixels)
    %   Instead of rows/cols we use a general resolution variable. This
    %   is because the eye model can only take equal rows and columns
    %   and the rendered image is always square.
    resolution;

    % fov - Field of view of the render in degrees
    %   This value is calculated from the retina distance and the
    %   retina size. This is only a close approximation since the
    %   retina is very slightly curved.
    fov;

    % accommodation - Diopters of accommodation for 550 nm light
    %   We change the properties of the lens to match the desired
    %   accommodation. For example, if we set this to 5 diopters, 550
    %   nm rays from 0.2 meters will be in focus on the retina.
    accommodation;

    % eccentricity - [Currently not implemented!] Horizontal and vertical
    %   angles on the retina corresponding to the center of the rendered
    %   image. Positive angles are to the right/up (from the eye's point of
    %   view) and negative angles are to the left/down. For example, an
    %   image with [0 0] eccentricity is centered on the center of the
    %   retina. An image with [30 0] eccentricity is centered 30 degrees to
    %   the right of the center of the retina.
    eccentricity;

    % pupilDiameter - Diameter of the pupil (mm)
    pupilDiameter;

    %retinaRadius - The curvature of the retina in mm
    %   If one imagines the retina as asection of a sphere, this radius
    %   value determines the distance from the edge of the sphere to
    %   its center. We will not change this most of the time, but
    %   sometimes it is helpful to make the retina very flat in order
    %   to measure certain properties of the eye.
    retinaRadius;

    % retinaDistance - Distance between the back lens and the retina
    %   We will not change this most of the time, but sometimes it is
    %   helpful to move the retina back and forth, like a camera
    %   sensor, to see things affects like chromatic aberration.
    retinaDistance;

    % numRays - Number of rays to shoot per pixel.
    %   This determines the quality of the render and affects the time
    %   spent rendering. This should be a factor of 2. Low quality is
    %   typically 64 or256 rays, high quality is typically 2048 or 4096
    %   rays.
    numRays;

    % numBounces - Number of bounces before ray terminates
    %   This also determines how accurately light is modeled in the
    %   rendering. The amount needed is highly scene dependent.
    %   Typically set to 1 for simple, diffuse scenes. A high value
    %   would be 4-8 for scenes with lots of reflections, caustics, or
    %   glassy materials.
    numBounces;

    % numCABands - Number of wavelength samples to take when modeling CA
    %   We shoot extra rays of different wavelengths in order to model
    %   chromatic aberration through the lens system. When debugging, this
    %   can be set to 0 but for the final render it should be something
    %   like 8 or 16. (e.g. If you set it to 8, then we will shoot rays for
    %   wavelengths of linspace(400, 700, 8).);
    numCABands;

    % eyePos - Position of the eye within the scene in [x y z] format
    %   [x y z]
    eyePos;

    % eyeTo - Point where the eye is looking at
    %   [x y z], the difference between eyeTo and eyePos is the
    %   direction vector that the optical axis is aligned with.
    eyeTo;

    % eyeUp - Up vector used when building the LookAt transform
    %   [x y z], this is typically [0 0 1] but it depends on how the
    %   eye is oriented. For example, if this was [0 0 -1] the eye
    %   would be "upside down." Some values are not valid, for example
    %   if the eye is looking down the z-axis (eyePos = [0 0 0], eyeTo
    %   = [0 0 1]) then the up vector cannot be [0 0 1].
    eyeUp;
    
    % diffractionEnabled Toggle diffraction. When it diffraction simulation
    %   is enabled, PBRT will use HURB to simulate the effect of
    %   diffraction. May cause slow-downs. 
    diffractionEnabled;
    
    %DEBUGMODE Toggle debug mode.
    %   For debug mode we switch to a perspective camera with the same
    %   FOV as the eye. This can be potentially faster and easier to
    %   render than going through the eye.
    debugMode;
    
    %RECIPE Structure that holds all other instructions needed for the
    %renderer
    % (PBRT) to render the scene.This includes things like the
    % WorldBegin/WorldEnd block, the PixelFilter, the Integrator, etc.
    % Ideally, the sceneEye user will not need to access the recipe very
    % often.
    recipe;
    
    % LENSFILE - Path to the .dat file that describes the lens system
    %   This file includes descriptions of the thickness, curvature, 
    %   and diameter of the various components in the eye. This is usually
    %   written out automatically in the "write" function for sceneEye
    %   (which is called during sceneEye.render). However, you can also
    %   attach a custom file. 
    lensFile;
    
    % LENSDENSITY - Lens pigment density. Equivalent to lens density in the
    % Lens class for ISETBio (non-3D calculations).
    lensDensity;
    
end

properties (Dependent)
    % width - Width of imaged retina (mm)
    %   Depends on fov, retinaDistance, and rows/cols
    width;

    % height - Height of imaged retina (mm)
    %   Depends on fov, retinaDistance, and rows/cols
    height;

    % sampleSize - Samples spacing, e.g. width/xRes and height/yRes.
    %   We assume square samples. This is not always accurate at large
    %   fov's.
    sampleSize;
    
    % angularSupport - location of each pixel in degrees. This should be
    % accurate even at wide-angles. May not be accurate if you use a crop
    % window though!
    angularSupport;

end

properties(GetAccess=public, SetAccess=private)

    % pbrtFile - Path to the original .pbrt file this scene is based on
    %   Depends on the pbrt file used to create the scene. Should not
    %   be changed.
    pbrtFile;

    % workingDir - Directory to store temp files needed for rendering
    %   We make a copy of the scene into the working directory, and
    %   then output new PBRT files into this directory. We also save
    %   the raw rendered data (xxx.dat) in this folder.
    workingDir;
    
    %SCENEUNITS Some scenes are in units of meters, some in units of millimeters.
    %   We keep of track of this here so we can pass the correct parameter
    %   to PBRT.
    sceneUnits;
      

end

properties(GetAccess=public, SetAccess=public, Hidden=true)
    
    %DISTANCE2CHORD This is used in intermediate calculations and is an
    %   important variable when we are doing calculations at wide angles.
    %   This is equivalent to the distance "a" shown in the diagram in
    %   get.width
    distance2chord;

end

properties (Constant)
    % wave - In PBRT we samples from 400 to 700 nm in 31 intervals
    wave = linspace(400, 700, 31); % nm

end

methods
    % Constructor
    function obj = sceneEye(pbrtFile, varargin)
        % Initialize the sceneEye class
        %
        % Reads a PBRT file and fills out the information needed
        % for the sceneEye object. That object will be rendered
        % using the PBRT methods (docker image).

        p = inputParser;
        p.KeepUnmatched = true;
        
        % pbrtFile: Either a pbrt file or just a scene name
        p.addRequired('pbrtFile', @ischar);
        p.addParameter('name', 'scene-001', @ischar);
        p.addParameter('workingDirectory', '', @ischar);
        
        % Parse
        p.parse(pbrtFile, varargin{:});
        
        % Setup the pbrt scene and recipe
        [recipe, obj.sceneUnits, obj.workingDir, obj.pbrtFile]  = ...
            loadPbrtScene(pbrtFile, varargin);
        
        % Check to make sure this PBRT file has a realistic eye.
        % [Note: JNM - 5/14/19 the human eye has retinaDistance parameters,
        % realistic does not. Changing type to see if call for ChessSet
        % continues failing.]
        if(~strcmp(recipe.camera.subtype, 'realisticEye'))
            % recipe.camera = piCameraCreate('realisticEye');
            % recipe.camera = piCameraCreate('realistic');
            recipe.camera = piCameraCreate('humaneye');
        end

        % Set properties
        obj.name = p.Results.name;
        obj.modelName = 'Navarro'; % Default
        obj.resolution = recipe.film.xresolution.value;
        obj.retinaDistance = recipe.camera.retinaDistance.value;
        obj.pupilDiameter = recipe.camera.pupilDiameter.value;

        obj.retinaDistance = recipe.camera.retinaDistance.value;
        obj.retinaRadius = recipe.camera.retinaRadius.value;

        retinaSemiDiam = recipe.camera.retinaSemiDiam.value;
        obj.fov = 2 * atand(retinaSemiDiam / obj.retinaDistance);

        % There's no variable for accommodation but we can infer it
        % from the name of the lens. We assume the naming conventions
        % is "%s_%f.dat" This is not foolproof, so maybe we can think
        % of a more robust way to do this in the future?
        obj.lensFile = recipe.camera.lensfile.value;
        if(strcmp(obj.lensFile, ''))
            obj.accommodation = [];
        else
            % Use regular expressions to find any floats within the string
            value = regexp(obj.lensFile, '(\d+, )*\d+(\.\d*)?', 'match');
            obj.accommodation = str2double(value{1});
        end

        obj.numRays = recipe.sampler.pixelsamples.value;

        % These two are often empty, so let's do checks here. However, 
        % I should find a more permanant solution to cases like these.
        % (See note above).
        % Maybe in piGetRenderRecipe we should put in the default
        % values if any of these rendering options are missing
        % (e.g. if Renderer is missing, put in Renderer 'sampler'.)
        if(isfield(recipe.integrator, 'maxdepth'))
            obj.numBounces = recipe.integrator.maxdepth.value;
        else
            obj.numBounces = 1;
        end
        if(isfield(recipe.renderer, 'nWaveBands'))
            obj.numCABands = recipe.renderer.nWaveBands.value;
        else
            obj.numCABands = 0;
        end

        if(~isempty(recipe.lookAt))
            obj.eyePos = recipe.lookAt.from;
            obj.eyeTo = recipe.lookAt.to;
            obj.eyeUp = recipe.lookAt.up;
        end

        obj.recipe = recipe;
        
        % Default settings.
        obj.debugMode = false;
        obj.diffractionEnabled = false;
        obj.eccentricity = [0 0];
        obj.lensDensity = 1.0;
        
    end

    %% Get methods for dependent variables
    % In PBRT, the image height is equivalent to the size of the chord on
    % the back of the spherical retina. We have to do the complex
    % calculation below in order to find an image size that would give us
    % the desired FOV. We want to measure the FOV from the back of the
    % lens. 
    %
    % I will attempt to illustrate this in ascii:
    %{
                       ooo OOO OOO ooo
                   oOO                 OOo
               oOO                         OOo
            oOO                               OOo
          oOO                                   OOo
        oOO                                     | OOo
       oOO                                      |  OOo
      oOO       back of lens                    |   OOo
     oOO             |                          |     OOo
     oOO             |---d----|-------a --------|--b--OOo
     oOO-------------*--------x-----------------------OOo  retina
     oOO                      |                ||     OOo
     oOO                      |                ||    OOo
      oOO               center of sphere     x ||   OOo
       oOO                                     ||  OOo
        oOO                                    || OOo
          oOO                                  -OOo
            oO                                OOo
               oOO                         OOo
                   oOO                 OOo
                       ooo OOO OOO ooo
    %}
    % Calculations:
    % retinaDistance = d + a + b
    % a^2 + x^2 = r^2
    % x/(d+a) = tand(fov/2) = k
    %
    % x = sqrt(r^2-a^2)
    % sqrt(r^2-a^2)/(d+a) - k = 0 
    % We can solve for a using an fzero solve.

    function val = get.distance2chord(obj)
        % Not entirely accurate but lets treat the origin point for the FOV
        % calculate as the beack of the lens
        if(obj.retinaRadius > obj.retinaDistance)
            error('Retina radius is larger than retina distance.')
        end

        myfun = @(a, k, d, r) sqrt(r^2-a.^2)./(d+a) - k;  % parameterized function
        k = tand(obj.fov/2);
        d = obj.retinaDistance - obj.retinaRadius;
        r = obj.retinaRadius;

        fun = @(a) myfun(a, k, d, r);    % function of x alone
        a = fzero(fun, [d obj.retinaRadius]);

        if(isnan(a))
            error('Search for a image width to match FOV failed. Initial guess is probably not close...')
        end

        val = a+d;
    end

    function val = get.width(obj)
        % Rendered image is alway square.
        val = 2 * tand(obj.fov / 2) * (obj.distance2chord);
    end

    function val = get.height(obj)
        % Rendered image is alway square.
        val = obj.width;
    end

    function val = get.sampleSize(obj)
        val = obj.width / obj.resolution;
    end

    function val = get.angularSupport(obj)
       % We have to be careful with this calculation.
       % Conver the chord distances to accurate angles.
       ss = obj.sampleSize;
       chordSpatialSamples = (1:obj.resolution).*ss - obj.width/2;
       val = atand(chordSpatialSamples/obj.distance2chord);
    end

    %% Set methods for dependent variables

    % Does this go here? MATLAB doesn't like this setup, but I would like
    % retinaDistance and retinaRadius to be both dependent (changes with
    % modelName), but also set-able by the user. What's the best way to do
    % this?

    % When we set the eye model, we need to change the retina distance and
    % radius.
    function set.modelName(obj, val)
        switch val
            case {'Navarro', 'navarro'}
                obj.modelName = 'Navarro';
                obj.retinaDistance = 16.32;
                obj.retinaRadius = 12;
            case {'LeGrand', 'legrand', 'le grand'}
                obj.modelName = 'LeGrand';
                obj.retinaDistance = 16.6;
                obj.retinaRadius = 13.4;
            case {'Arizona', 'arizona'}
                obj.modelName = 'Arizona';
                obj.retinaDistance = 16.713;
                obj.retinaRadius = 13.4;
            case {'Custom', 'custom'}
                % Use Navarro as a default.
                % Do we need to be able to change this later?
                obj.modelName = 'Custom';
                if(isempty(obj.retinaDistance))
                    obj.retinaDistance = [];
                end
                if(isempty(obj.retinaRadius))
                    obj.retinaRadius = [];
                end
        end
    end

    % When the user toggles debugMode, make sure the camera type is
    % correct.
    function set.debugMode(obj, val)
        obj.debugMode = val;
        if(val)
            obj.modelName = 'none';
            % The camera will be changed to perspective in write(), so we
            % do nothing here. 
        elseif(~val && strcmp(obj.modelName, 'none'))
            % Put the navarro eye back in if there's not already a model.
            obj.modelName = 'Navarro';
            obj.recipe.camera = piCameraCreate('humaneye');
        end
    end

    % I want to put in this warning, but again MATLAB doesn't really like
    % this!
    function set.accommodation(obj, val)
        obj.accommodation = val;
        if(strcmp(obj.modelName, 'Gullstrand'))
            warning(['Gullstrand eye has no accommodation modeling.', ...
                'Setting accommodation will do nothing.']);
        end
    end

    function set.lensFile(obj, val)
        
        % On creation, the lensFile is left empty
        % There should be a better way to do this right? I don't think I'm
        % doing this right. Maybe we need a seperate set function like we
        % do with render recipes?
        if(~isempty(val))
            % Make sure it's a valid file
            [p, ~, e] = fileparts(val);
            if(isempty(p))
                error('Lens file needs to be a full file path.');
            elseif(~strcmp(e, '.dat'))
                error('Lens file needs to be a .dat file.');
            end

            obj.modelName = 'Custom';
            obj.lensFile = val;
        else
            obj.lensFile = '';
        end
    end
    
    % I don't think these are necessary. 
%{
    function set.width(obj, val)
        obj.fov = 2 * atand((val / 2) / obj.retinaDistance);
    end

    function set.height(obj, val)
        obj.fov = 2 * atand((val / 2) / obj.retinaDistance);
    end

    function set.sampleSize(obj, val)
        obj.fov = 2 * atand((val * obj.resolution / 2) ...
            / obj.retinaDistance);
    end
%}
end

methods (Access=public)
    [oi, terminalOutput, outputFile] = render(obj, varargin);
    
    % These are helper functions called within render() above. Splitting
    % them into their individual functions allows us to integrate them with
    % isetcloud tools.
    % (Should these be public?)
    [obj] = setOI(obj, ieObject, varargin)
    [objNew] = write(obj, varargin)
end

end
