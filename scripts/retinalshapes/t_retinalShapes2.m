%% t_retinalShapes2
%
% Integrate
% See also
%   t_cones*

%% Initialize and clear
ieInit;

%% Build a simple scene and oi (retinal image) for computing

% First the scene
s = sceneCreate('rings rays');
s = sceneSet(s, 'fov', 1);

% Then the oi
oi = oiCreate;
oi = oiCompute(oi, s);

%% Build a default cone mosaic and compute isomerizatoins

% Get the default set of cone mosaic parameters.
cmParams = cMosaicParams;

% Generate off-axis mosaic
cmParams.eccentricityDegs = [0 0];   % The size of the cones out here is about the same as the mouse?
cmParams.sizeDegs = [1.5 1.5]*2;
cmParams.micronsPerDegree = oiGet(oi,'distance per degree','um');
cm = cMosaic(cmParams);

cm.visualize;

%% Printout the JSON file for TG's new implementation

% The file should include coneRFPositionsMicrons and a Z-value
% We should include the Z-computation Matlab script that TG wrote in
% ISETBio.
%
% We will need to integrate for different amounts of area on the
% rendered OI to account for cone aperture size.  Or will we?  Also,
% what assumption do we need to make about the orientation of the
% aperture?  SCE to think about.
%

% These are the positions.  
pos = cm.coneRFpositionsMicrons;
pos = pos*1e-3;   % Positions in millimeters?  Or meters?


%% Load the computed z-values

% or use the spherical map or other tools to create the z-values.


%% Compute isomerizations for each eye position.
[noiseFree, noisy] = cm.compute(oi);

vParams = cm.visualize('params');
vParams.activation = noisy;
vParams.activationColorMap = gray(512);
vParams.verticalActivationColorBar = true;
vParams.activationRange = [0 max(noiseFree(:))];
cm.visualize(vParams);

%% END
