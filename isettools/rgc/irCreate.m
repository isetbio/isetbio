function obj = irCreate(inputObj, varargin)
% Generate an ir object.
%
%  obj = irCreate(inputObj, params)
% 
%           params:
%             'name', name
%             'eyeSide',{'left','right'},
%             'eyeRadius',eyeRadius
%             'eyeAngle',eyeAngle
%             
% Inputs: 
%  inputObj:  a bipolar object or osDisplayRGBo object
%  name:      Name for this instance
%  model:     Computation type for all the rgc mosaics
%       'linear' - Basic linear filtering as in typical center-surround
%       'LNP'    - linear-nonlinear-Poisson, see Pillow paper as below; only contains post-spike filter
%       'GLM'    - coupled generalized linear model, see Pillow paper, includes coupling filters  
%       'Phys'   - pulls a set of parameters measured in physiology by
%                           the Chichilnisky lab.
%
% Outputs: 
%  rgc object: of the specified model type
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
% See also:  ir.m, rgcMosaic.m
%
% JRG 9/2015 Copyright ISETBIO Team
% JRG 7/2016 updated

%% Parse inputs
p = inputParser;
p.addRequired('inputObj');
addParameter(p,'name','inner retina',@ischar);
addParameter(p,'species','unknown',@ischar);

% In the future, we will read these from the input object
addParameter(p,'eyeSide',    'left', @ischar);
addParameter(p,'eyeRadius',   4,     @isnumeric);
addParameter(p,'eyeAngle',    0,     @isnumeric);  % X-axis is 0, positive Y is 90

p.parse(inputObj,varargin{:});

%% Create the object
obj = ir(inputObj, p.Results);


