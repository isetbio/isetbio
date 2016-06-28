function obj = irCreate(os, varargin)
%% irCreate: generate an @rgcLinear, @rgcLNP or @rgcGLM object.
%
%  obj = irCreate(outersegment,'name',name,'eyeSide',{'left','right'},'eyeRadius',eyeRadius,'eyeAngle',eyeAngle)
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
%   os  = osCreate('identity');
%   innerRetina = irCreate(os,'GLM','name','myRGC'); 
%   innerRetina = irCreate(os,'name','EJ',...
%             'eyeSide','left','eyeRadius',12,'eyeAngle',90));
% 
% See also:  t_rgc.m
%
% JRG 9/2015 Copyright ISETBIO Team

%% Parse inputs
p = inputParser;
p.addRequired('os');
addParameter(p,'name','inner retina',@ischar);
addParameter(p,'species','unknown',@ischar);

% In the future, we will read these from the os object, not here.  JRG is
% adding these parameters
addParameter(p,'eyeSide',    'left', @ischar);
addParameter(p,'eyeRadius',   4,     @isnumeric);
addParameter(p,'eyeAngle',    0,     @isnumeric);  % X-axis is 0, positive Y is 90

p.parse(os,varargin{:});

% Maybe delete?? BW/JRG
switch class(os)
    case{'osIdentity'}
        % For fast testing.  Creates a random image
        if isempty(osGet(os,'rgbData'))
            os = osSet(os, 'rgbData', rand(64));
        end
        
        if isempty(osGet(os,'patch size'))
            os = osSet(os, 'patch size', 180);
        end
        
        if isempty(osGet(os,'time step'))
            os = osSet(os,'time step',.01);
        end

    otherwise
        % Don't worry, carry on
end

%% Create the object
obj = ir(os, p.Results);


