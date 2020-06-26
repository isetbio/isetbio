function sSupport = oiSpatialSupport(oi, units)
% Calculate the spatial positions of the optical image sample points
%
% Syntax:
%   sSupport = oiSpatialSupport(oi, [units])
%
% Description:
%    Determine the spatial support for the optical image. The positions are
%    specified in x and y positions measured on the surface of the image
%    sensor. The units (meters, mm, um) can be specified.
%
%    There are examples contained in the source code below. To access, type
%    'edit oiSpatialSupport.m' into the Command Window.
%
% Inputs:
%    oi       - Struct. The optical image structure.
%    units    - (Optional) String. A string specifying units that support
%               is returned in. Default 'meters'. See oiGet(oi, 'spatial
%               resolution') for info on available units.
%
% Outputs:
%    sSupport - Struct with fields x and y, each specifing the spatial
%               support of the corresponding dimension of the oi.
%
% Optional key/value pairs:
%    None.
%

% History:
%	 xx/xx/03       Copyright ImagEval Consultants, LLC, 2003.
%    03/06/18  jnm  Formatting
%    06/24/19  JNM  Documentation pass.

% Examples:
%{
   scene = sceneCreate;
   oi = oiCreate;
   oi = oiCompute(oi, scene);
   sSupportmm = oiSpatialSupport(oi, 'millimeters')
   sSupportum = oiSpatialSupport(oi, 'microns')
   sSupportm = oiSpatialSupport(oi)
%}

if notDefined('units'), units = 'meters'; end

sr = oiGet(oi, 'spatialResolution', units);
nRows = oiGet(oi, 'rows');
nCols = oiGet(oi, 'cols');

if isempty(nRows) || isempty(nCols), error('No optical image data.'); end

sSupport.y = linspace(-nRows * sr(1) / 2 + sr(1) / 2, ...
    nRows * sr(1) / 2 - sr(1) / 2, nRows);
sSupport.x = linspace(-nCols * sr(2) / 2 + sr(2) / 2, ...
    nCols * sr(2) / 2 - sr(2) / 2, nCols);

end