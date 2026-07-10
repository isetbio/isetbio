%% t_eyeDiffraction.m
% SkipFile
% Multi-render human-eye diffraction comparison is useful interactively, but
% currently slow/unstable for routine smoke tests (CPU/remote dependent).
%
% We recommend you go through t_eyeIntro.m before running
% this tutorial.
%
% This tutorial renders a retinal image of "slanted edge." We can then use
% this slanted edge to estimate the modulation transfer function of the
% optical system.
%
% We also show how the color fringing along the edge of the bar due to
% chromatic aberration. 
%
% Depends on: 
%   pbrt2ISET, ISETBIO, Docker, ISET
%
% TL ISETBIO Team, 2017

%% Check ISETBIO and initialize

ieInit;
if ~piDockerExists, piDockerConfig; end

%% Set up the slanted edge scene

thisSE = sceneEye('slanted edge');

thisSE.set('rays per pixel',32);
thisSE.set('render type',{'radiance','depth'});

from = [0 0 -100];
thisSE.set('from',from);
thisSE.get('lookat')

thisSE.set('use pinhole',true);

scene = thisSE.piWRS('name','pinhole');

% piAssetGeometry(thisSE.recipe);

%% Runs, but units are not set right.  Too blurry.
% Maybe accommodation.

thisSE.set('use optics',true);

thisSE.set('fov',1);                % About 3 deg on a side
thisSE.set('spatial samples',256);  % Number of OI sample points
thisSE.set('rays per pixel',256);
% tmp = 1/thisSE.get('object distance','m');
thisSE.set('accommodation',1);  
thisSE.set('lens density',0);       % Yellow is harder to see.

thisSE.set('diffraction',true);
thisSE.set('pupil diameter',4);
thisSE.set('film diagonal',10);

oi = thisSE.piWRS('name','navarro');

%%
oi = oiSet(oi,'name','4mm-diffractionOn');
oiWindow(oi);
oiPlot(oi,'illuminance hline',[128 128]);

%% Diffraction should not matter

thisSE.set('diffraction',false);
thisSE.set('rays per pixel',512);
thisSE.set('pupil diameter',1);

oi = thisSE.piWRS('name','navarro');
oiPlot(oi,'illuminance hline',[128 128]);
title('4 mm off')

%% Diffraction should matter

thisSE.set('diffraction',true);
oi = thisSE.piWRS('name','navarro-diffraction');

oiPlot(oi,'illuminance hline',[128 128]);
title('1 mm on')
%% Diffraction should matter.

% Make a direct comparison
thisSE.set('diffraction',false);
oi = oiSet(oi,'name','1mm-diffractionOff');
oiWindow(oi);
thisSE.summary;

oiPlot(oi,'illuminance hline',[128 128]);
set(gca,'xlim',[-30 30],'xtick',(-30:10:30));
title('1 mm off')
%%  Maybe we should be smoothing the curve at the edge?

thisSE.set('rays per pixel',4096);
thisSE.set('pupil diameter',0.5);

thisSE.set('diffraction',true);
oi = oiSet(oi,'name','Halfmm-diffractionOn');
oiWindow(oi);

oiPlot(oi,'illuminance hline',[128 128]);
set(gca,'xlim',[-30 30],'xtick',(-30:10:30));
title('Half mm on')
%%

thisSE.set('diffraction',false);
oi = oiSet(oi,'name','Halfmm-diffractionOff');
oiWindow(oi);

oiPlot(oi,'illuminance hline',[128 128]);
set(gca,'xlim',[-30 30],'xtick',(-30:10:30));
title('Half mm off')
%% END
