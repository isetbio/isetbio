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
        [scene, display] = sceneHorwitzHassWhiteNoise(params);
        oi  = oiCreate('wvf human');
        sensor = sensorHorwitzHassShortWhiteNoise(params, scene, oi, display);
        identityOS = osCreate('identity');
        sceneRGB = sceneHorwitzHassWhiteNoiseRGB(params);               
        identityOS = osSet(identityOS, 'rgbData', sceneRGB);
        obj = rgcGLM(sensor, identityOS, 'right', 3.0, 180);
        
    elseif nargin == 6
        if strcmpi(varargin{1},'linear');
            obj = rgcLinear(varargin{2:6});
        elseif strcmpi(varargin{1},'lnp');
            obj = rgcLNP(varargin{2:6});
        else strcmpi(varargin{1},'glm');
            obj = rgcGLM(varargin{2:6});
        end
        
    else
        warning(sprintf('\nrgcCreate: create an isetbio @rgc object.\nrgcCreate(type, sensor, outersegment, eyeSide, patchRadius, patchAngle)'));
    end
    