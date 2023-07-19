%% Introduction to the cone mosaic (cMosaic) object.
%
% Description:
%    Create a cone mosaic object and compute cone isomerizations.
%    No eye movements.
%
%    Visualize the results with a plot.
%
% See also
%   t_cones*

%% Initialize and clear
ieInit;

%% Build a simple scene and oi (retinal image) for computing

% First the scene
s = sceneCreate('rings rays');
s = sceneSet(s, 'fov', 1);

% Then the oi
oi = oiCreate('human');
oi = oiCompute(oi, s);

%% Build a default cone mosaic and compute isomerizatoins

% Get the default set of cone mosaic parameters.
cmParams = cMosaicParams;

cmParams.eccentricityDegs = [0 0];
cmParams.sizeDegs = [1.5 1.5];
cmParams.micronsPerDegree = oiGet(oi,'distance per degree','um');
cm = cMosaic(cmParams);

cm.visualize;

%% Compute isomerizations for each eye position.
[noiseFree, noisy] = cm.compute(oi);

vParams = cm.visualize('params');
vParams.activation = noisy;
vParams.activationColorMap = gray(512);
vParams.verticalActivationColorBar = true;
vParams.activationRange = [0 max(noiseFree(:))];
cm.visualize(vParams);

%% END
