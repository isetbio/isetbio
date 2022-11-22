classdef sceneEye < hiddenHandle
% Create a sceneEye object
%
% Syntax:
%   thisSE = sceneEye(sceneName,varargin);
%
% Inputs
%  sceneName - 
% 
% Optional key/val
%    'eye model' - 'navarro','legrand','arizona'
%     
% Output
%    sceneEye - Modified scene eye object
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
%    The sceneEye includes the PBRT rendering recipe (thisR) in one of
%    its slots. The camera struct is called 'human eye' and it
%    contains the slots that are necessary to specify the human eye
%    model.  These parameter slots differ from the standard camera
%    model (e.g., 'realistic', 'pinhole', or 'omni').  The slots in
%    'human eye' include retinal curvature, the index of refraction
%    of the components of the eye, and so forth.  These are fixed, and
%    to change accommodation we change the lens model using
%    recipe.set('accommodation') and related functions. 
%
%    Dependencies: ISET3D
%
% See Also:
%    t_eyeNavarro, setNavarroAccommodation, recipe

% Examples:
%{
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
        % Initialize the sceneEye class.  Default human eye model is
        % navarro.
        %
        %   thisSE = sceneEye;
        %   thisSE = sceneEye('bedroom');
        %   thisSE = sceneEye('bathroom','human eye','legrand');
        %
        % Reads a PBRT file and fills out the information needed
        % for the sceneEye object. That object will be rendered
        % using the PBRT methods (docker image).

        if ~exist('pbrtFile','var'), pbrtFile = ''; end
        varargin = ieParamFormat(varargin);

        p = inputParser;
        p.addRequired('pbrtfile',@(x)(isempty(x) || ischar(x)));
        p.addParameter('eyemodel','navarro',@(x)ismember(x,{'navarro','legrand','arizona'}));
        p.parse(pbrtFile,varargin{:});

        % Setup the pbrt scene recipe
        if isempty(pbrtFile),  obj.recipe = recipe;
        else,                  obj.recipe = piRecipeDefault('scene name',pbrtFile);
        end
        
        % Create the camera model
        obj.set('camera',piCameraCreate('humaneye','eye model',p.Results.eyemodel));
        
        % At this point the camera is created.  The recipe should have an
        % output dir, so we can create the default lens file
        switch (p.Results.eyemodel)
            case 'navarro'
                navarroWrite(obj.recipe,0);
            case 'arizona'
                arizonaWrite(obj.recipe,0);
            case 'legrand'
                legrandWrite(obj.recipe);
        end        

        % Assign this object the basename of the input file
        obj.set('name',obj.get('input basename'));
                        
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
