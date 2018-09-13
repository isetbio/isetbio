% Script for testing the oiPlot routine
%
% Description:
%    Exercises various plot options available for the optical image
%    structure, via routine oiPlot. There are more options available than
%    shown here, but this gives you the idea.
%
% See Also:
%    oiPlot, oiGet.
%

% History
%    12/30/17  dhb  Added comments, changed name from scripts/s_oiPlot ->
%                   tutorials/t_optics/t_oiPlot
%    09/10/18  jnm  Formatting

%% Initialize ISETBio and the OI structure
ieInit;

scene = sceneCreate; 
scene = sceneSet(scene, 'fov', 4);
oi = oiCreate('diffraction limited'); 
oi = oiCompute(oi, scene);

%% Irradiadiance along a vertical line
% Plotted as a function of position and wavelength.
%
% Here the vector roiLocs gives x, y position. Only the x position is used
% when plotting a vertical line through the image.
roiLocs = [20 20];
[uData, g] = oiPlot(oi, 'irradiance vline', roiLocs);

%% Irradiadiance along a horizontal line
% Plotted as a function of position and wavelength.
%
% Here only the y position of roiLocs is used.
[uData, g] = oiPlot(oi, 'irradiance hline', roiLocs);

% Now along the middle row of the image. 
middleRow = round(oiGet(oi, 'rows') / 2);
uData = oiPlot(oi, 'irradiance hline', [1, middleRow]);

%% Retinal illuminance along a horiztonal line
[uData, g] = oiPlot(oi, 'illuminance hline', roiLocs);

% Transform of illuminance along the middle row
uData = oiPlot(oi, 'illuminance fft hline', [1, middleRow]);

%% Contrast along a horizontal line
uData = oiPlot(oi, 'contrast hline', [1, middleRow]);

%% Retinal Irradiance
% Image of retinal irradiance
uData = oiPlot(oi, 'irradiance image with grid', [], 40);

% Retinal irradiance at 500 nm, with a grid at 40 um spacing
uData = oiPlot(oi, 'irradiance image wave', [], 500, 40);

% Transform of retinal irradiance at 450 nm
uData = oiPlot(oi, 'irradiance fft', [], 450);

%% Transform of retinal illuminance
uData = oiPlot(oi, 'illuminance fft');

%%  Average spectrum over a roi
% This involves user interaction, because the roi is not passed explicitly,
% so it is commented out.
%
% Select and execute.
%{
    uData = oiPlot(oi, 'irradiance energy roi');
%}

%% PSF at 550 nm, plotted as a function of um
uData = oiPlot(oi, 'psf 550', 'um');

%% OTF at 550 nm, plotted as a funciton of um
uData = oiPlot(oi, 'otf 550', 'um');

%% Line spread function versus wavelength
uData = oiPlot(oi, 'ls wavelength');
