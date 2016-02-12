function obj = irCreate(os, varargin)
%% irCreate: generate an @rgcLinear, @rgcLNP or @rgcGLM object.
%
%  obj = irCreate(outersegment,'name',name,'model',{linear,LNP,GLM,Subunit,phys})
% 
% Inputs: 
%  os:     a outer segment structure
%  name:   Name for this instance
%  model:   Computation type for all the rgc mosaics
%       'linear' - Basic linear filtering as in typical center-surround
%       'LNP'    - linear-nonlinear-Poisson, see Pillow paper as below; only contains post-spike filter
%       'GLM'    - coupled generalized linear model, see Pillow paper, includes coupling filters  
%       'Subunit'- cones act as individual subunits
%       'Phys'   - pulls a set of parameters measured in physiology by
%                           the Chichilnisky lab.
%
% Outputs: 
%  rgc object: There are several types divided according to the model
%              implementation.
% 
% The coupled-GLM model is described in Pillow, Shlens, Paninski, Sher,
% Litke, Chichilnisky & Simoncelli, Nature (2008).
% 
% The LNP and GLM models here are based on
% <http://pillowlab.princeton.edu/code_GLM.html code by Pillow> 
% under the GNU General Public License.
%  
% Example:
%   os  = osCreate('linear');
%   rgc = irCreate(os,'GLM','name','myRGC','type','LNP'); 
%   rgc = irCreate(os,'type','GLM', 'name','EJ');
% 
% See also:  The initialize method for the @rgcLinear, @rgcLNP or @rgcGLM
%            subclasses define the specific implementations.
%
% TODO:  Deal with eye parameter issues (see below).
%
% JRG 9/2015 Copyright ISETBIO Team

%% Parse inputs
p = inputParser;
p.addRequired('os');
addParameter(p,'name','inner retina',@ischar);
addParameter(p,'model','linear',@ischar);
addParameter(p,'species','unknown',@ischar);

% In the future, we will read these from the os object, not here.  JRG is
% adding these parameters
addParameter(p,'eyeSide',    'left', @ischar);
addParameter(p,'eyeRadius',   4,     @isnumeric);
addParameter(p,'eyeAngle',    0,     @isnumeric);  % X-axis is 0, positive Y is 90

p.parse(os,varargin{:});

%% Create the object
switch ieParamFormat(p.Results.model)
    case {'linear','rgclinear'}
        obj = rgcLinear(os, p.Results);
    case {'pool','rgcpool'}
        obj = rgcPool(os, p.Results);
    case {'lnp','rgclnp'}
        obj = rgcLNP(os, p.Results);
    case {'glm','rgcglm'}
        obj = rgcGLM(os, p.Results);
    case {'subunit','rgcsubunit'}        
        obj = rgcSubunit(os, p.Results);
    case{'phys','rgcphys'}
        % NEED TO SWITCH ORDER HERE!
        obj = rgcPhys(outersegment, sensor, 'eyeSide', eyeSide, 'eyeRadius', eyeRadius, 'eyeAngle', eyeAngle);
    otherwise
        error('Unrecognized model type, please choose from ''linear'', ''LNP'', ''GLM'', etc.');
        
end

