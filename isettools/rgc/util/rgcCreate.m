function obj = rgcCreate(type, varargin)
% function obj = rgcCreate(type, sensor, outersegment, eyeSide, patchRadius, patchAngle)
% rgcCreate: generate an @rgcLinear, @rgcLNP or @rgcGLM object.
% 
%   rgc = rgcCreate(model type, params)
%         where params = [...
%                 scene, ...
%                 sensor, ...
%                 os, ...
%                 eyeSide, ...
%                 eyeRadius, ...
%                 eyeAngle];
% 
% Inputs: 
%       model type: 
%            'linear' - 
%            'LNP'    - linear-nonlinear-Poisson, see Pillow paper as below; only contains post-spike filter
%            'GLM'    - coupled generalized linear model, see Pillow paper, includes coupling filters  
%            'Subunit'- cones act as individual subunits
%            'Phys'   - pulls a set of parameters measured in physiology by
%                           the Chichilnisky lab.
% 
%       scene:  a scene structure
%       sensor: a sensor structure
%       os:     a outer segment structure
%
%    Optional, but highly recommended:
%       eyeSide: 'left' or 'right', which eye the retinal patch is from
%       patchRadius: radius of retinal patch in microns
%       patchAngle: polar angle of retinal patch
%     [These inputs determine the size of spatial receptive fields, and are
%       necessary to accurately model physiological responses.]
% 
% Outputs: 
%   rgc object
% 
% The coupled-GLM model is described in Pillow, Shlens, Paninski, Sher,
% Litke, Chichilnisky & Simoncelli, Nature (2008).
% 
% The LNP and GLM models here are based on code by Pillow available at 
%       http://pillowlab.princeton.edu/code_GLM.html
% under the GNU General Public License.
%  
% See the initialize method for the @rgcLinear, @rgcLNP or @rgcGLM 
% subclasses for more details of the specific implementations.
%
% Example:
% % default values assumed for eyeSide, eyeRadius, eyeAngle.
%   rgc2 = rgcCreate('GLM', scene, sensor, os); 
%
%   eyeAngle = 180; % degrees
%   eyeRadius = 3; % mm
%   eyeSide = 'right';
%   rgc1 = rgcCreate('GLM', scene, absorptions, os, eyeSide, eyeRadius, eyeAngle);
%
% 
% JRG 9/2015

%%

narginchk(0, Inf);
p = inputParser; p.CaseSensitive = false; p.FunctionName = mfilename;

addRequired( p, 'type');
addParameter(p,'scene',        'sc1');
addParameter(p,'sensor',       's1');
addParameter(p,'outersegment', 'os1');
addParameter(p,'eyeSide',      'side', @ischar);
addParameter(p,'eyeRadius',   8,     @isnumeric);
addParameter(p,'eyeAngle',    0,     @isnumeric);

% Parse and put results into structure p.
p.parse(type, varargin{:}); 
scene        = p.Results.scene;
sensor       = p.Results.sensor;
outersegment = p.Results.outersegment;

% Eye parameters
eyeSide   = p.Results.eyeSide;
eyeRadius = p.Results.eyeRadius;
eyeAngle  = p.Results.eyeAngle;

% Set key-value pairs.
switch lower(p.Results.type)
    case {'linear','rgclinear'}
        obj = rgcLinear(scene, sensor, outersegment, eyeSide, eyeRadius, eyeAngle);
    case {'lnp','rgclnp'}
        obj = rgcLNP(scene, sensor, outersegment, eyeSide, eyeRadius, eyeAngle);
    case {'glm','rgcglm'}
        obj = rgcGLM(scene, sensor, outersegment, eyeSide, eyeRadius, eyeAngle);
    case {'subunit','rgcsubunit'}
        obj = rgcSubunit(scene, sensor, outersegment, eyeSide, eyeRadius, eyeAngle);
    case{'phys','rgcphys'}
        obj = rgcPhys(scene, sensor, outersegment, eyeSide, eyeRadius, eyeAngle);
    otherwise
        obj = rgcLinear(scene, sensor, outersegment, eyeSide, eyeRadius, eyeAngle);
  
end

