function sSupport = sceneSpatialSupport(scene, units)
% Calculate the spatial positions of the scene sample locations
%
% Syntax:
%	sSupport = sceneSpatialSupport(scene, [units])
%
% Description:
%    Determine the spatial positions of the sample positions (spatial
%    support) of a scene. These values correspond to position the scene.
%    The default spatial units are returned as part of a structure in x and
%    y positions in meters.
%
%    N.B. The source contains executable examples of usage, which can be
%    accessed by typing 'edit sceneSpatialSupport.m' in the command window.
%
% Inputs:
%    scene    - The scene structure
%    units    - (Optional) The measurement units. Default is 'meters'.
%
% Outputs:
%    sSupport - The calculated scene spatial positions
%
% Optional key/value pairs:
%    None.
%
% See Also:
%    oiSpatialSupport
%

% History:
%    xx/xx/03       Copyright ImagEval Consultants, LLC, 2003.
%    12/20/17  jnm  Formatting & fix example
%    01/25/18  jnm  Formatting update to match Wiki.

% Examples:
%{
    scene = sceneCreate;
    sSupportmm = sceneSpatialSupport(scene, 'millimeters');
    sSupportum = sceneSpatialSupport(scene, 'microns');
    sSupportm  = sceneSpatialSupport(scene);
%}

if notDefined('units'), units = 'meters'; end

sr = sceneGet(scene, 'spatialResolution', units);
nRows = sceneGet(scene, 'rows');
nCols = sceneGet(scene, 'cols');

sSupport.y = linspace(-nRows * sr(1) / 2 + sr(1) / 2, ...
    nRows * sr(1) / 2 - sr(1) / 2, nRows);
sSupport.x = linspace(-nCols * sr(2) / 2 + sr(2) / 2, ...
    nCols * sr(2) / 2 - sr(2) / 2, nCols);
          
end