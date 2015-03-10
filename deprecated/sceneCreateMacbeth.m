function scene = sceneCreateMacbeth(surface,lightSource,scene)
% Create a hyperspectral scene of the Macbeth chart.
%   
%   scene = sceneCreateMacbeth(surface,lightSource,[scene])
%
% The surface reflectances and light source are specified as arguments.
% The color signal is computed (in photons); these values are then
% attached to the scene structure.
%
% Used by the sceneWindow callbacks to create the Macbeth images.
%
% Copyright ImagEval Consultants, LLC, 2005.

warning('Deprecated. Use sceneCreate instead');

if notDefined('lightSource.spectrum')
    error('Bad light source description.');
elseif notDefined('surface.spectrum')
    error('Bad surface description.');
elseif ~isequal(lightSource.spectrum.wave(:),surface.spectrum.wave(:))
    error('Mis-match between light source and object spectral fields')
end

% This is particularly unclear code because of the way the Macbeth surfaces
% are stored on disk.  That should change.  I believe what is happening is
% we take the light source in photons and we make it the same size and
% shape as the surface data.  Then we multiply point by point.
iPhotons = illuminantGet(lightSource,'photons');
[surface,r,c] = RGB2XWFormat(surface.data);
sPhotons = surface*diag(iPhotons);
sPhotons = XW2RGBFormat(sPhotons,r,c);

% We compute the product of the surface reflectance and illuminant photons
% here
% scene = sceneSet(scene,'cphotons',surface.data .* photons);
scene = sceneSet(scene,'photons',sPhotons);

% Store the light source
scene = sceneSet(scene,'illuminant',lightSource);

end