classdef sceneEye < hiddenHandle
% Create a sceneEye object
%
% Syntax:
%   myScene = sceneEye();
%
% Description:
%    sceneEye is represents the information needed to construct a new PBRT
%    file that we render to estimate the retinal image (spectral
%    irradiance) of the human eye.
%
%    The sceneEye class is analogous to the "scene" structure in ISETBIO
%    (and ISET). The class supports similar methods. Unlike the older
%    ISET/ISETBIO "scene" struct, sceneEye is implemented as a MATLAB class
%    with its own methods.
%
%    The sceneEye includes the PBRT rendering recipe (thisR) in one of its
%    slots. The camera struct is called 'realisticEye' and it contains the
%    slots that are necessary to specify the human eye model.  These
%    parameter slots differ from the standard camera model (e.g.,
%    'realistic', 'pinhole', or 'omni'.  The slots in 'realisticEye'
%    include retinal curvature, the index of refraction of the components
%    of the eye, and so forth.
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
%    Dependencies: ISET3D
%
% See Also:
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
    
    % USEPINHOLE - a debugging mode
    %   For debug mode we switch to a pinhole camera with the same
    %   FOV as the eye. This can be potentially faster and easier to
    %   render than going through the eye.
    usePinhole = false;
    
    %RECIPE Structure that holds all other instructions needed for the
    %renderer
    % (PBRT) to render the scene.This includes things like the
    % WorldBegin/WorldEnd block, the PixelFilter, the Integrator, etc.
    % Ideally, the sceneEye user will not need to access the recipe very
    % often.
    recipe;
    
    % LENSDENSITY - Lens pigment density. Equivalent to lens density in the
    % Lens class for ISETBio (non-3D calculations).
    lensDensity = 1;
    
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
      
end

properties(GetAccess=public, SetAccess=public, Hidden=true)
    

end

properties (Constant)
    % wave - In PBRT we samples from 400 to 700 nm in 31 intervals
    % wave = linspace(400, 700, 31); % nm

end

methods
    % Constructor
    function obj = sceneEye(pbrtFile, varargin)
        % Initialize the sceneEye class
        %
        % Reads a PBRT file and fills out the information needed
        % for the sceneEye object. That object will be rendered
        % using the PBRT methods (docker image).

        if ~exist('pbrtFile','var'), pbrtFile = ''; end
         
        % Setup the pbrt scene recipe
        if isempty(pbrtFile),  obj.recipe = recipe;
        else,                  obj.recipe = piRecipeDefault('scene name',pbrtFile);
        end
        
        % Assign this object the basename of the input file
        obj.set('name',obj.get('input basename'));
        
        % Make sure the recipe specifies realistic eye.  That camera has
        % the parameters needed to model the human.  Note:  realisticEye
        % differs from realistic.
        % disp('Setting Navarro model default');
        obj.set('camera',piCameraCreate('humaneye','lens file','navarro.dat'));
        obj.modelName = 'navarro'; % Default
                
    end

    %% Get methods for dependent variables
    
    % In PBRT, the image height is equivalent to the size of the chord on
    % the back of the spherical retina. We have to do the calculation below
    %  to find an image size that would give us the desired FOV. We
    % want to measure the FOV from the back of the lens.
    %
    % I will attempt to illustrate this in ascii - and there is a
    % PowerPoint in the iset3d.wiki/images directory that does this, too.
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
    % a^2 + x^2 = r^2   (radius of the sphere)
    % x/(d+a) = tand(fov/2) = k
    %
    % x = sqrt(r^2-a^2)
    % sqrt(r^2-a^2)/(d+a) - k = 0 
    % We can solve for a using an fzero solve.

    % These are pretty much constant unless the model changes them.
    % r = retinaRadius
    % d + a + b = retinaDistance
    %
    % This one should determine the field of view, partitioning a and b.
    %
    % x = retinaSemiDiam

    
    %{
    function val = get.distance2chord(obj)
        % Not entirely accurate but lets treat the origin point for the FOV
        % calculate as the back of the lens
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
    %}

    % Get function for sceneEye. Returns derived parameters of the
    % sceneEye that require some computation
    function val = get(obj,param, varargin)
        val = eyeGet(obj, param, varargin{:});
    end
    
    % Sets parameters of sceneEye. Mostly these are passed through to the
    % rendering recipe.
    function obj = set(obj,param,varargin)
        obj = eyeSet(obj,param,varargin{:});
    end
    
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
