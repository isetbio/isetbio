function scene = sceneInit(varargin)
% First stage in building up a scene
%
% Syntax
%   scene = sceneInit(...);
%
% Description
%   These are the first few lines of sceneCreate.  Sometimes we want them
%   available for testing other functions, such as sceneHarmonic and the
%   like.
%
% Input
%   None
%
% Return
%   scene struct
%
% Optional key/val pairs
%   bit depth - Default is 32, which is single precision
% 
% BW, ISETBIO Team, 2017
%
% See also
%    sceneCreate

%{
  scene = sceneInit;
  scene = sceneInit('bit depth',64);
  scene = sceneInit('bit depth',32);
%}
p = inputParser;
p.addParameter('bitdepth',32,@(x)(x == 32 || x == 64));  % Single precision by default

varargin = ieParamFormat(varargin);
p.parse(varargin{:});
bitdepth = p.Results.bitdepth;

% Identify the object type
scene.type = 'scene';

% Set precision
scene = sceneSet(scene,'bit depth',bitdepth);

end