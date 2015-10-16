function obj = rgcCreate(varargin)
% rgcCreate: generate an @rgcLinear, @rgcLNP or @rgcGLM object.
% 
% Inputs: none.
% Outputs: the rgc object.
% 
% See the initialize method for the @rgcLinear, @rgcLNP or @rgcGLM 
% subclasses for more details of the specific implementations.
%
% 9/2015 JRG


    if nargin == 0
        params.image_size = 64; params.meanLuminance = 100;
        params.nsteps = 30; params.fov = 0.8;
        fprintf(sprintf('\nGenerating scene, oi and display:\n'));
        [scene, display] = sceneHorwitzHassWhiteNoise(params);
        oi  = oiCreate('wvf human');
        sensor = sensorHorwitzHassShortWhiteNoise(params, scene, oi, display);
        identityOS = osCreate('identity');
        sceneRGB = sceneHorwitzHassWhiteNoiseRGB(params);               
        identityOS = osSet(identityOS, 'rgbData', sceneRGB);
        obj = rgcGLM(scene, sensor, identityOS, 'right', 3.0, 180);
        
    elseif nargin == 1
        params.image_size = 64; params.meanLuminance = 100;
        params.nsteps = 30; params.fov = 0.8;
        fprintf(sprintf('\nGenerating scene, oi and display:\n'));
        [scene, display] = sceneHorwitzHassWhiteNoise(params);
        oi  = oiCreate('wvf human');
        sensor = sensorHorwitzHassShortWhiteNoise(params, scene, oi, display);
        identityOS = osCreate('identity');
        sceneRGB = sceneHorwitzHassWhiteNoiseRGB(params);               
        identityOS = osSet(identityOS, 'rgbData', sceneRGB);
        % obj = rgcGLM(scene, sensor, identityOS, 'right', 3.0, 180);
        
        if strcmpi(varargin{1},'linear');
            obj = rgcLinear(scene, sensor, identityOS, 'right', 3.0, 180);
        elseif strcmpi(varargin{1},'lnp');
            obj = rgcLNP(scene, sensor, identityOS, 'right', 3.0, 180);
        elseif strcmpi(varargin{1},'glm');
            obj = rgcGLM(scene, sensor, identityOS, 'right', 3.0, 180);
        else strcmpi(varargin{1},'subunit');
            obj = rgcSubunit(scene, sensor, identityOS, 'right', 3.0, 180);
        end
        
    elseif nargin == 7
        if strcmpi(varargin{1},'linear')|strcmpi(varargin{1},'rgclinear');
            obj = rgcLinear(varargin{2:7});
        elseif strcmpi(varargin{1},'lnp')|strcmpi(varargin{1},'rgclnp');
            obj = rgcLNP(varargin{2:7});
        elseif strcmpi(varargin{1},'glm')|strcmpi(varargin{1},'rgcglm');
            obj = rgcGLM(varargin{2:7});
        else strcmpi(varargin{1},'subunit')|strcmpi(varargin{1},'rgcsubunit');
            obj = rgcSubunit(varargin{2:7});
        end
        
    else
        warning(sprintf('\nrgcCreate: create an isetbio @rgc object.\nrgcCreate(type, sensor, outersegment, eyeSide, patchRadius, patchAngle)'));
    end
    