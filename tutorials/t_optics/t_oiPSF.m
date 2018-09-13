% PSF plotting and manipulation from the oi and optics
%
% Description:
%    Plot psfs at 500 nm for the diffraction limited, wvf human, and
%    marimont-wandell models. They differ quite a bit in one sense. But the
%    spatial support of the wvf and MW models are somewhat similar.
%
% Notes:
%    * [Note: XXX - Programming note: I am not happy that the spatial
%      sampling of the PSF is so arbitrary and poorly controlled. We need
%      to figure out when the frequency/spatial sampling gets set and do a
%      better job of controlling it.]
%
% See Also:
%    t_oiPlot

% History:
%    xx/xx/18       Copyright, ISETBIO Team, 2018
%    09/10/18  jnm  Formatting

%% Initialization & Diffraction test
ieInit

% Test for diffraction case
scene = sceneCreate;

%% Diffraction limited optics
oi = oiCreate('diffraction');
oi = oiCompute(oi, scene);

oiPlot(oi, 'psf550');
psf = oiGet(oi, 'optics psf data', 500);
support = oiGet(oi, 'optics psf support', 'um');
mesh(support{1}, support{2}, psf)
xlabel('um');
ylabel('um');
set(gca, 'xlim', [-25 25], 'ylim', [-25 25])

%% Wavefront model
oi = oiCreate('wvf human');
oi = oiCompute(oi, scene);

oiPlot(oi, 'psf550');
psf = oiGet(oi, 'optics psf data', 500);
support = oiGet(oi, 'optics psf support', 'um');
mesh(support{1}, support{2}, psf)
xlabel('um');
ylabel('um');
set(gca, 'xlim', [-25 25], 'ylim', [-25 25])

%% Marimont and Wandell human model
oi = oiCreate('human');
oi = oiCompute(oi, scene);

oiPlot(oi, 'psf550');
psf = oiGet(oi, 'optics psf data', 500);
support = oiGet(oi, 'optics psf support', 'um');
mesh(support{1}, support{2}, psf)
xlabel('um');
ylabel('um');
set(gca, 'xlim', [-25 25], 'ylim', [-25 25])
