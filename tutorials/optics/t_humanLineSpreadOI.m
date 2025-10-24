% Demonstrate retinal irradiance of line stimuli at varying wavelengths
%
% Description:
%    This function is intended to illustrate the retinal irradiance of a
%    line stimuli at different wavelengths. 
%

% History:
%    XX/XX/15       Copyright ISETBIO Team, 2015
%    11/26/18  JNM  Formatting

%% Initialize
ieInit;

%% Create a line scene, human optics, and a human sensor
% This is a broad band stimulus, with a spectral power distribution of
% daylight, 6500 K.  We set the field of view to one degree.
lineS = sceneCreate('line d65', [128, 256]);
lineS = sceneSet(lineS, 'h fov', 1);
% sceneWindow(lineS);

% In the past, we used the Marimont and Wandell estimates from
% the mid-90s. 
% oi = oiCreate('human mw');

% These optics are the estimated human optics from Thibos' mean.
oi = oiCreate('human wvf');

%%  Compute and display the broad band

oi = oiCompute(oi,lineS,'pad value','mean','crop',true);

oiWindow(oi);

roi = [];
wList = [450, 550, 650];  % nm
gSpacing = 40;            % microns
figHdl = ieFigure([],'wide');
tiledlayout(1,3)
for ww = 1:length(wList)
    thisWave = wList(ww);
    nexttile;
    oiPlot(oi, 'irradiance image wave grid', roi, thisWave, gSpacing, figHdl);    

    % Make the grid lines more visible
    ax = get(figHdl,'CurrentAxes');
    ax.GridLineWidth = 3;
    ax.GridColor = [0.8 0.8 0.8];
end

%% END