%% s_sceneDataExtractionAndPlotting
%
% It is possible to interactively select regions of interest in the scene
% window and plot the scene energy, photons, or reflectance.  The
% reflectance is derived by dividing the scene photons and the illuminant
% (in photons).
%
% In general when we open ISET scene data we either read an illuminant or we
% assume a default (D65).   We set the illuminant level to be consistent
% with a peak reflectance term.  This method is a little imprecise right
% now because we don't specify the wavelength of the peak reflectance.  We
% will fix this some day soon.
%
% Copyright ImagEval Consultants, LLC, 2010

%% Create a scene and plot reflectance
scene = sceneCreate('macbethd65');
vcAddAndSelectObject(scene);
sceneWindow;  % The user selects the region to plot interactively

% Here are the luminance data from a line
rows = round(sceneGet(scene,'rows')/2);
[uData, h] = plotScene(scene,'luminance hline',[1,rows]);

%% The uData  structure  

% This structure contains the data in the graph
uData

% This structure is also attached to the figure
get(h,'userdata')

%% Scenes store information about the illuminant
plotScene(scene,'illuminant energy roi')

%% You can plot the energy 
rect = [51    35    10    11];        % Yellow Macbeth patch
roiLocs = ieRoi2Locs(rect);  % xy locations in scene
plotScene(scene,'radiance energy roi',roiLocs);

%% Mean quanta (photons) at the corresponding location
plotScene(scene,'radiance photons roi',roiLocs);

%% You can also estimate the reflectance
plotScene(scene,'reflectance',roiLocs);

%% If you would just like to read the data with no plot, you can use
radiance = vcGetROIData(scene,roiLocs,'photons');
radiance = mean(radiance);
wave = sceneGet(scene,'wave');
vcNewGraphWin;
plot(wave,radiance)
ylabel('Radiance (q/s/nm/m^2/sr)');

%% Or for energy
radiance = vcGetROIData(scene,roiLocs,'photons');
radiance = Quanta2Energy(wave,radiance);
radiance = mean(radiance);
vcNewGraphWin;
plot(wave,radiance)
ylabel('Radiance (watts/s/m^2/nm/sr)');

%% End
