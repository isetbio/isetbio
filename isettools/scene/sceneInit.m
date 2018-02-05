function scene = sceneInit(varargin)
% First stage in building up a scene
%
% Syntax
%   scene = sceneInit([varargin]);
%
% Description
%   These are the first few lines of sceneCreate.  Sometimes we want them
%   available for testing other functions, such as sceneHarmonic and the
%   like.
%
%    There are examples in the code. Type 'edit sceneInit' into the Command
%    Window to access.
%
% Inputs:
%    None required.
%
% Outputs:
%    scene    - The scene structure
%
% Optional key/val pairs:
%    bitDepth - The bit depth. Default is 32, which is single precision
%
% See Also:
%    scenecreate

% History:
%    xx/xx/17  BW   ISETBIO Team, 2017
%    01/25/18  jnm  Formatting

% Examples:
%{
  scene = sceneInit;
  scene = sceneInit('bit depth', 64);
  scene = sceneInit('bit depth', 32);
%}
p = inputParser;
% Single precision bit depth by default
p.addParameter('bitdepth', 32, @(x)(x == 32 || x == 64));

varargin = ieParamFormat(varargin);
p.parse(varargin{:});
bitdepth = p.Results.bitdepth;

% Identify the object type
scene.type = 'scene';

% Set precision
scene = sceneSet(scene, 'bit depth', bitdepth);

end