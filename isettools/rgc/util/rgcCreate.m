function obj = rgcCreate(varargin)
% rgcCreate: generate an @rgcLinear, @rgcLNP or @rgcGLM object.
%
%  obj = rgcCreate(type, sensor, outersegment, eyeSide, patchRadius, patchAngle)
% 
% The arguments can be sent in as name/value pairs, or the parser
% allows us to send in the parameters as elements of a structure,
% say
%
%  REDO THESE ARGUMENTS!!!!
%  params = [...
%           sensor, ... (TODO:  row,col of cones and spacing)
%           os, ...
%           eyeSide, ...
%           eyeRadius, ...
%           eyeAngle ...
%   ];
% We are still figuring out what is exactly needed here
%
% Inputs: 
%  model type: 
%       'linear' - 
%       'LNP'    - linear-nonlinear-Poisson, see Pillow paper as below; only contains post-spike filter
%       'GLM'    - coupled generalized linear model, see Pillow paper, includes coupling filters  
%       'Subunit'- cones act as individual subunits
%       'Phys'   - pulls a set of parameters measured in physiology by
%                           the Chichilnisky lab.
%       sensor: a sensor structure
%       os:     a outer segment structure
%
%  Optional, but highly recommended:
%       eyeSide: 'left' or 'right', which eye the retinal patch is from
%       patchRadius: radius of retinal patch in microns
%       patchAngle: polar angle of retinal patch
%       [These inputs determine the size of spatial receptive fields, and are
%       necessary to accurately model physiological responses.]
% 
% Outputs: 
%   rgc object
% 
% The coupled-GLM model is described in Pillow, Shlens, Paninski, Sher,
% Litke, Chichilnisky & Simoncelli, Nature (2008).
% 
% The LNP and GLM models here are based on
% <http://pillowlab.princeton.edu/code_GLM.html code by Pillow> 
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

addParameter(p,'name','inner retina',@ischar);
addParameter(p,'model','linear',@ischar);
addParameter(p,'species','unknown',@ischar);

addParameter(p,'row',     [], @isnumeric);
addParameter(p,'col',     [], @isnumeric);
addParameter(p,'spacing', [], @isnumeric);
addParameter(p,'timing', [], @isnumeric);

addParameter(p,'eyeSide',      'side', @ischar);
addParameter(p,'eyeRadius',   8,     @isnumeric);
addParameter(p,'eyeAngle',    0,     @isnumeric);

p.parse(varargin{:});

% Set key-value pairs
switch ieParamFormat(p.Results.model)
    case {'linear','rgclinear'}
        obj = rgcLinear(p.Results);
        % obj = rgcLinear(outersegment, sensor, 'eyeSide', eyeSide, 'eyeRadius', eyeRadius, 'eyeAngle', eyeAngle);
    case {'lnp','rgclnp'}
        obj = rgcLNP(outersegment, sensor, 'eyeSide', eyeSide, 'eyeRadius', eyeRadius, 'eyeAngle', eyeAngle);
    case {'glm','rgcglm'}
        obj = rgcGLM(p.Results);
        % obj = rgcGLM(outersegment, sensor, 'eyeSide', eyeSide, 'eyeRadius', eyeRadius, 'eyeAngle', eyeAngle);
    case {'subunit','rgcsubunit'}
        obj = rgcSubunit(outersegment, sensor, 'eyeSide', eyeSide, 'eyeRadius', eyeRadius, 'eyeAngle', eyeAngle);
    case{'phys','rgcphys'}
        % NEED TO SWITCH ORDER HERE!
        obj = rgcPhys(outersegment, sensor, 'eyeSide', eyeSide, 'eyeRadius', eyeRadius, 'eyeAngle', eyeAngle);
    otherwise
        obj = rgcLinear(outersegment, sensor, 'eyeSide', eyeSide, 'eyeRadius', eyeRadius, 'eyeAngle', eyeAngle);
  
end

